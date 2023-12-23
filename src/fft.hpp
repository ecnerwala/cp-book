#pragma once

#include <vector>
#include <cmath>
#include <utility>
#include <algorithm>

#include "modnum.hpp"

/**
 * Author: Andrew He
 * Source: http://neerc.ifmo.ru/trains/toulouse/2017/fft2.pdf
 * Papers about accuracy: http://www.daemonology.net/papers/fft.pdf, http://www.cs.berkeley.edu/~fateman/papers/fftvsothers.pdf
 * For integers rounding works if $(|a| + |b|)\max(a, b) < \mathtt{\sim} 10^9$, or in theory maybe $10^6$.
 */

namespace ecnerwala {
namespace fft {

using std::swap;
using std::vector;
using std::min;
using std::max;

template<class T> int sz(T&& arg) { using std::size; return int(size(std::forward<T>(arg))); }
inline int nextPow2(int s) { return 1 << (s > 1 ? 32 - __builtin_clz(s-1) : 0); }

// Complex
template <typename dbl> struct cplx { /// start-hash
	dbl x, y;
	cplx(dbl x_ = 0, dbl y_ = 0) : x(x_), y(y_) { }
	friend cplx operator+(cplx a, cplx b) { return cplx(a.x + b.x, a.y + b.y); }
	friend cplx operator-(cplx a, cplx b) { return cplx(a.x - b.x, a.y - b.y); }
	friend cplx operator*(cplx a, cplx b) { return cplx(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
	friend cplx conj(cplx a) { return cplx(a.x, -a.y); }
	friend cplx inv(cplx a) { dbl n = (a.x*a.x+a.y*a.y); return cplx(a.x/n,-a.y/n); }
};

// getRoot implementations
template <typename num> struct getRoot {
	static num f(int k) = delete;
};
template <typename dbl> struct getRoot<cplx<dbl>> {
	static cplx<dbl> f(int k) {
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
		dbl a=2*M_PI/k;
		return cplx<dbl>(cos(a),sin(a));
	}
};
template <int MOD> struct primitive_root {
	static const int value;
};
template <> struct primitive_root<998244353> {
	static const int value = 3;
};
template <int MOD> struct getRoot<modnum<MOD>> {
	static modnum<MOD> f(int k) {
		assert((MOD-1)%k == 0);
		return pow(modnum<MOD>(primitive_root<MOD>::value), (MOD-1)/k);
	}
};

template <typename num> class fft {
	static vector<int> rev;
	static vector<num> rt;

public:
	static void init(int n);
	template <typename Iterator> static void go(Iterator begin, int n);

	static vector<num> scratch_a;
	static vector<num> scratch_b;
};

template <typename num> vector<int> fft<num>::rev;
template <typename num> vector<num> fft<num>::rt;
template <typename num> vector<num> fft<num>::scratch_a;
template <typename num> vector<num> fft<num>::scratch_b;

template <typename num> void fft<num>::init(int n) {
	if (n <= sz(rt)) return;
	rev.resize(n);
	for (int i = 0; i < n; i++) {
		rev[i] = (rev[i>>1] | ((i&1)*n)) >> 1;
	}
	rt.reserve(n);
	while (sz(rt) < 2 && sz(rt) < n) rt.push_back(num(1));
	for (int k = sz(rt); k < n; k *= 2) {
		rt.resize(2*k);
		num z = getRoot<num>::f(2*k);
		for (int i = k/2; i < k; i++) {
			rt[2*i] = rt[i], rt[2*i+1] = rt[i]*z;
		}
	}
}

template <typename num> template <typename Iterator> void fft<num>::go(Iterator begin, int n) {
	init(n);
	int s = __builtin_ctz(sz(rev)/n);
	for (int i = 0; i < n; i++) {
		if (i < (rev[i]>>s)) {
			swap(*(begin+i), *(begin+(rev[i]>>s)));
		}
	}
	for (int k = 1; k < n; k *= 2) {
		for (int i = 0; i < n; i += 2 * k) {
			Iterator it1 = begin + i, it2 = it1 + k;
			for (int j = 0; j < k; j++, ++it1, ++it2) {
				num t = rt[j+k] * *it2;
				*it2 = *it1 - t;
				*it1 = *it1 + t;
			}
		}
	}
}

template <typename num> struct fft_multiplier {
	template <typename IterA, typename IterB, typename IterOut>
	static void multiply(IterA ia, int sza, IterB ib, int szb, IterOut io) {
		vector<num>& fa = fft<num>::scratch_a;
		vector<num>& fb = fft<num>::scratch_b;

		if (sza == 0 || szb == 0) return;
		int s = sza + szb - 1;
		int n = nextPow2(s);
		if (sz(fa) < n) fa.resize(n);
		if (sz(fb) < n) fb.resize(n);
		fft<num>::init(n);

		bool did_cut = false;
		if (sza > 1 && szb > 1 && n == 2 * (s - 1)) {
			// we have exactly 1 wraparound, so let's just handle it explicitly to save a factor of 2
			// only do it if sza < s and szb < s so we don't have to wrap the inputs
			did_cut = true;
			n /= 2;
		}
		copy(ia, ia+sza, fa.begin());
		fill(fa.begin()+sza, fa.begin()+n, num(0));
		copy(ib, ib+szb, fb.begin());
		fill(fb.begin()+szb, fb.begin()+n, num(0));
		// used if did_cut
		num v_init; if (did_cut) { v_init = fa[0] * fb[0]; }
		fft<num>::go(fa.begin(), n);
		fft<num>::go(fb.begin(), n);
		num d = inv(num(n));
		for (int i = 0; i < n; i++) fa[i] = fa[i] * fb[i] * d;
		reverse(fa.begin()+1, fa.begin()+n);
		fft<num>::go(fa.begin(), n);
		if (did_cut) {
			fa[s-1] = std::exchange(fa[0], v_init) - v_init;
		}
		copy(fa.begin(), fa.begin()+s, io);
	}

	template <typename IterA, typename IterOut>
	static void square(IterA ia, int sza, IterOut io) {
		multiply<IterA, IterA, IterOut>(ia, sza, ia, sza, io);
	}
};

template <typename num>
struct fft_inverser {
	template <typename IterA, typename IterOut>
	static void inverse(IterA ia, int sza, IterOut io) {
		vector<num>& fa = fft<num>::scratch_a;
		vector<num>& fb = fft<num>::scratch_b;

		if (sza == 0) return;
		int s = nextPow2(sza) * 2;
		fft<num>::init(s);
		if (sz(fa) < s) fa.resize(s);
		if (sz(fb) < s) fb.resize(s);
		fb[0] = inv(*ia);
		for (int n = 1; n < sza; ) {
			fill(fb.begin() + n, fb.begin() + 4 * n, num(0));
			n *= 2;
			copy(ia, ia+min(n,sza), fa.begin());
			fill(fa.begin()+min(n,sza), fa.begin()+2*n, 0);
			fft<num>::go(fb.begin(), 2*n);
			fft<num>::go(fa.begin(), 2*n);
			num d = inv(num(2*n));
			for (int i = 0; i < 2*n; i++) fb[i] = fb[i] * (2 - fa[i] * fb[i]) * d;
			reverse(fb.begin()+1, fb.begin()+2*n);
			fft<num>::go(fb.begin(), 2*n);
		}
		copy(fb.begin(), fb.begin()+sza, io);
	}
};

template <typename dbl>
struct fft_double_multiplier {
	template <typename IterA, typename IterB, typename IterOut>
	static void multiply(IterA ia, int sza, IterB ib, int szb, IterOut io) {
		vector<cplx<dbl>>& fa = fft<cplx<dbl>>::scratch_a;
		vector<cplx<dbl>>& fb = fft<cplx<dbl>>::scratch_b;

		if (sza == 0 || szb == 0) return;
		int s = sza + szb - 1;
		int n = nextPow2(s);
		fft<cplx<dbl>>::init(n);
		if (sz(fa) < n) fa.resize(n);
		if (sz(fb) < n) fb.resize(n);

		fill(fa.begin(), fa.begin() + n, 0);
		{ auto it = ia; for (int i = 0; i < sza; ++i, ++it) fa[i].x = *it; }
		{ auto it = ib; for (int i = 0; i < szb; ++i, ++it) fa[i].y = *it; }
		fft<cplx<dbl>>::go(fa.begin(), n);
		for (auto& x : fa) x = x * x;
		for (int i = 0; i < n; ++i) fb[i] = fa[(n-i)&(n-1)] - conj(fa[i]);
		fft<cplx<dbl>>::go(fb.begin(), n);
		{ auto it = io; for (int i = 0; i < s; ++i, ++it) *it = fb[i].y / (4*n); }
	}

	template <typename IterA, typename IterOut>
	static void square(IterA ia, int sza, IterOut io) {
		multiply<IterA, IterA, IterOut>(ia, sza, ia, sza, io);
	}
};

template <typename mnum>
struct fft_mod_multiplier {
	template <typename IterA, typename IterB, typename IterOut>
	static void multiply(IterA ia, int sza, IterB ib, int szb, IterOut io) {
		using cnum = cplx<double>;
		vector<cnum>& fa = fft<cnum>::scratch_a;
		vector<cnum>& fb = fft<cnum>::scratch_b;

		if (sza == 0 || szb == 0) return;
		int s = sza + szb - 1;
		int n = nextPow2(s);
		fft<cnum>::init(n);
		if (sz(fa) < n) fa.resize(n);
		if (sz(fb) < n) fb.resize(n);

		{ auto it = ia; for (int i = 0; i < sza; ++i, ++it) fa[i] = cnum(int(*it) & ((1<<15)-1), int(*it) >> 15); }
		fill(fa.begin()+sza, fa.begin() + n, 0);
		{ auto it = ib; for (int i = 0; i < szb; ++i, ++it) fb[i] = cnum(int(*it) & ((1<<15)-1), int(*it) >> 15); }
		fill(fb.begin()+szb, fb.begin() + n, 0);

		fft<cnum>::go(fa.begin(), n);
		fft<cnum>::go(fb.begin(), n);
		double r0 = 0.5 / n; // 1/2n
		for (int i = 0; i <= n/2; i++) {
			int j = (n-i)&(n-1);
			cnum g0 = (fb[i] + conj(fb[j])) * r0;
			cnum g1 = (fb[i] - conj(fb[j])) * r0;
			swap(g1.x, g1.y); g1.y *= -1;
			if (j != i) {
				swap(fa[j], fa[i]);
				fb[j] = fa[j] * g1;
				fa[j] = fa[j] * g0;
			}
			fb[i] = fa[i] * conj(g1);
			fa[i] = fa[i] * conj(g0);
		}
		fft<cnum>::go(fa.begin(), n);
		fft<cnum>::go(fb.begin(), n);
		using ll = long long;
		const ll m = mnum::MOD;
		auto it = io;
		for (int i = 0; i < s; ++i, ++it) {
			*it = mnum((ll(fa[i].x+0.5)
						+ (ll(fa[i].y+0.5) % m << 15)
						+ (ll(fb[i].x+0.5) % m << 15)
						+ (ll(fb[i].y+0.5) % m << 30)) % m);
		}
	}

	template <typename IterA, typename IterOut>
	static void square(IterA ia, int sza, IterOut io) {
		multiply<IterA, IterA, IterOut>(ia, sza, ia, sza, io);
	}
};

template <class multiplier, typename num>
struct multiply_inverser {
	template <typename IterA, typename IterOut>
	static void inverse(IterA ia, int sza, IterOut io) {
		if (sza == 0) return;
		int s = nextPow2(sza);
		vector<num> b(s,num(0));
		vector<num> tmp(2*s);
		b[0] = inv(*ia);
		for (int n = 1; n < sza; ) {
			multiplier::square(b.begin(),n,tmp.begin());
			int nn = min(sza,2*n);
			multiplier::multiply(tmp.begin(),nn,ia,nn,tmp.begin());
			for (int i = n; i < nn; i++) b[i] = -tmp[i];
			n = nn;
		}
		copy(b.begin(), b.begin()+sza, io);
	}
};

template <class multiplier, typename T> vector<T> multiply(const vector<T>& a, const vector<T>& b) {
	if (sz(a) == 0 || sz(b) == 0) return {};
	vector<T> r(sz(a) + sz(b) - 1);
	multiplier::multiply(begin(a), sz(a), begin(b), sz(b), begin(r));
	return r;
}

template <typename T> vector<T> fft_multiply(const vector<T>& a, const vector<T>& b) {
	return multiply<fft_multiplier<T>, T>(a, b);
}
template <typename T> vector<T> fft_double_multiply(const vector<T>& a, const vector<T>& b) {
	return multiply<fft_double_multiplier<T>, T>(a, b);
}
template <typename T> vector<T> fft_mod_multiply(const vector<T>& a, const vector<T>& b) {
	return multiply<fft_mod_multiplier<T>, T>(a, b);
}

template <class multiplier, typename T> vector<T> square(const vector<T>& a) {
	if (sz(a) == 0) return {};
	vector<T> r(2 * sz(a) - 1);
	multiplier::square(begin(a), sz(a), begin(r));
	return r;
}
template <typename T> vector<T> fft_square(const vector<T>& a) {
	return square<fft_multiplier<T>, T>(a);
}
template <typename T> vector<T> fft_double_square(const vector<T>& a) {
	return square<fft_double_multiplier<T>, T>(a);
}
template <typename T> vector<T> fft_mod_square(const vector<T>& a) {
	return square<fft_mod_multiplier<T>, T>(a);
}

template <class inverser, typename T> vector<T> inverse(const vector<T>& a) {
	vector<T> r(sz(a));
	inverser::inverse(begin(a), sz(a), begin(r));
	return r;
}
template <typename T> vector<T> fft_inverse(const vector<T>& a) {
	return inverse<fft_inverser<T>, T>(a);
}
template <typename T> vector<T> fft_double_inverse(const vector<T>& a) {
	return inverse<multiply_inverser<fft_double_multiplier<T>, T>, T>(a);
}
template <typename T> vector<T> fft_mod_inverse(const vector<T>& a) {
	return inverse<multiply_inverser<fft_mod_multiplier<T>, T>, T>(a);
}
/* namespace fft */ }

// Power series; these are assumed to be the min of the length
template <typename T, typename multiplier, typename inverser>
struct power_series : public std::vector<T> {
	using std::vector<T>::vector;

	int ssize() const {
		return int(this->size());
	}
	int len() const {
		return ssize();
	}
	int degree() const {
		return len() - 1;
	}
	void extend(int sz) {
		assert(sz >= ssize());
		this->resize(sz);
	}
	void shrink(int sz) {
		assert(sz <= ssize());
		this->resize(sz);
	}
	// multiply by x^n
	void shift(int n = 1) {
		assert(n >= 0 && n <= ssize());
		std::rotate(this->begin(), this->end()-n, this->end());
		std::fill(this->begin(), this->begin()+n, T(0));
	}
	// divide by x^n and 0-pad
	void unshift(int n = 1) {
		assert(n >= 0 && n <= ssize());
		std::fill(this->begin(), this->begin()+n, T(0));
		std::rotate(this->begin(), this->begin()+n, this->end());
	}
	power_series& operator += (const power_series& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	friend power_series operator + (const power_series& a, const power_series& b) {
		power_series r(std::min(a.size(), b.size()));
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] + b[i];
		}
		return r;
	}
	power_series& operator -= (const power_series& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
	}
	friend power_series operator - (const power_series& a, const power_series& b) {
		power_series r(std::min(a.size(), b.size()));
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] - b[i];
		}
		return r;
	}

	power_series& operator *= (const T& n) {
		for (auto& v : *this) v *= n;
		return *this;
	}
	friend power_series operator * (const power_series& a, const T& n) {
		power_series r(a.size());
		for (int i = 0; i < a.len(); i++) {
			r[i] = a[i] * n;
		}
		return r;
	}
	friend power_series operator * (const T& n, const power_series& a) {
		power_series r(a.size());
		for (int i = 0; i < a.len(); i++) {
			r[i] = n * a[i];
		}
		return r;
	}

	friend power_series operator * (const power_series& a, const power_series& b) {
		if (sz(a) == 0 || sz(b) == 0) return {};
		power_series r(std::max(0, sz(a) + sz(b) - 1));
		multiplier::multiply(begin(a), sz(a), begin(b), sz(b), begin(r));
		r.resize(std::min(a.size(), b.size()));
		return r;
	}
	power_series& operator *= (const power_series& o) {
		return *this = (*this) * o;
	}
	friend power_series square(const power_series& a) {
		if (sz(a) == 0) return {};
		power_series r(sz(a) * 2 - 1);
		multiplier::square(begin(a), sz(a), begin(r));
		r.resize(a.size());
		return r;
	}

	friend power_series stretch(const power_series& a, int n) {
		power_series r(a.size());
		for (int i = 0; i*n < int(a.size()); i++) {
			r[i*n] = a[i];
		}
		return r;
	}
	friend power_series inverse(power_series a) {
		power_series r(sz(a));
		inverser::inverse(begin(a), sz(a), begin(r));
		return r;
	}
	friend power_series deriv_shift(power_series a) {
		for (int i = 0; i < a.len(); i++) {
			a[i] *= i;
		}
		return a;
	}
	friend power_series integ_shift(power_series a) {
		assert(a[0] == 0);
		T f = 1;
		for (int i = 1; i < int(a.size()); i++) {
			a[i] *= f;
			f *= i;
		}
		f = inv(f);
		for (int i = int(a.size()) - 1; i > 0; i--) {
			a[i] *= f;
			f *= i;
		}
		return a;
	}
	friend power_series deriv_shift_log(power_series a) {
		auto r = deriv_shift(a);
		return r * inverse(a);
	}
	friend power_series poly_log(power_series a) {
		assert(a[0] == 1);
		return integ_shift(deriv_shift_log(std::move(a)));
	}
	friend power_series poly_exp(power_series a) {
		assert(a.size() >= 1);
		assert(a[0] == 0);
		power_series r(1, T(1));
		while (r.size() < a.size()) {
			int n_sz = std::min(int(r.size()) * 2, int(a.size()));
			r.resize(n_sz);
			power_series v(a.begin(), a.begin() + n_sz);
			v -= poly_log(r);
			v[0] += 1;
			r *= v;
		}
		return r;
	}
	friend power_series poly_pow_monic(power_series a, int64_t k) {
		if (a.empty()) return a;
		assert(a.size() >= 1);
		assert(a[0] == 1);
		a = poly_log(a);
		a *= k;
		return poly_exp(a);
	}
	friend power_series poly_pow(power_series a, int64_t k) {
		int st = 0;
		while (st < a.len() && a[st] == 0) st++;
		if (st == a.len()) return a;

		power_series r(a.begin() + st, a.end());
		T leading_coeff = r[0];
		r *= inv(leading_coeff);
		r = poly_pow_monic(r, k);
		r *= pow(leading_coeff, k);
		r.insert(r.begin(), size_t(st), T(0));
		return r;
	}

	friend power_series to_newton_sums(const power_series& a, int deg) {
		auto r = log_deriv_shift(a);
		r[0] = deg;
		for (int i = 1; i < int(r.size()); i++) r[i] = -r[i];
		return r;
	}
	friend power_series from_newton_sums(power_series S, int deg) {
		assert(S[0] == int(deg));
		S[0] = 0;
		for (int i = 1; i < int(S.size()); i++) S[i] = -S[i];
		return poly_exp(integ_shift(std::move(S)));
	}

	// Calculates prod 1/(1-x^i)^{a[i]}
	friend power_series euler_transform(const power_series& a) {
		power_series r = deriv_shift(a);
		std::vector<bool> is_prime(a.size(), true);
		for (int p = 2; p < int(a.size()); p++) {
			if (!is_prime[p]) continue;
			for (int i = 1; i*p < int(a.size()); i++) {
				r[i*p] += r[i];
				is_prime[i*p] = false;
			}
		}
		return poly_exp(integ_shift(r));
	}
	friend power_series inverse_euler_transform(const power_series& a) {
		power_series r = deriv_shift(poly_log(a));
		std::vector<bool> is_prime(a.size(), true);
		for (int p = 2; p < int(a.size()); p++) {
			if (!is_prime[p]) continue;
			for (int i = (int(a.size())-1)/p; i >= 1; i--) {
				r[i*p] -= r[i];
				is_prime[i*p] = false;
			}
		}
		return integ_shift(r);
	}
};

template <typename num> using power_series_fft = power_series<num, fft::fft_multiplier<num>, fft::fft_inverser<num>>;
template <typename num, typename multiplier> using power_series_with_multiplier = power_series<num, multiplier, fft::multiply_inverser<multiplier, num>>;
template <typename num> using power_series_fft_mod = power_series_with_multiplier<num, fft::fft_mod_multiplier<num>>;
template <typename num> using power_series_fft_double = power_series_with_multiplier<num, fft::fft_double_multiplier<num>>;

// TODO: Use iterator traits to deduce value type?
template <typename base_iterator, typename value_type> struct add_into_iterator {
	base_iterator base;
	add_into_iterator() : base() {}
	add_into_iterator(base_iterator b) : base(b) {}
	add_into_iterator& operator * () { return *this; }
	add_into_iterator& operator ++ () { base.operator ++ (); return *this; }
	add_into_iterator& operator ++ (int) { auto temp = *this; operator ++ (); return temp; }
	auto operator = (value_type v) { base.operator * () += v; }
};

// TODO: Use iterator traits to deduce value type?
template <typename base_iterator, typename value_type> struct add_double_into_iterator {
	base_iterator base;
	add_double_into_iterator() : base() {}
	add_double_into_iterator(base_iterator b) : base(b) {}
	add_double_into_iterator& operator * () { return *this; }
	add_double_into_iterator& operator ++ () { base.operator ++ (); return *this; }
	add_double_into_iterator& operator ++ (int) { auto temp = *this; operator ++ (); return temp; }
	auto operator = (value_type v) { base.operator * () += 2 * v; }
};

template <typename num, typename multiplier> struct online_multiplier {
	int N; int i;
	std::vector<num> f, g;
	std::vector<num> res;

	// Computes the first 2N terms of the product
	online_multiplier(int N_) : N(N_), i(0), f(N), g(N), res(2*N+1, num(0)) {}

	num peek() {
		return res[i];
	}

	void push(num v_f, num v_g) {
		assert(i < N);
		f[i] = v_f;
		g[i] = v_g;
		if (i == 0) {
			res[i] += v_f * v_g;
		} else {
			res[i] += v_f * g[0];
			res[i] += f[0] * v_g;
			// TODO: We could do this second half more lazily, since it only affects res[i+1]...
			for (int p = 1; (i & (p-1)) == (p-1); p <<= 1) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				multiplier::multiply(
					// TODO: We can cache FFT([f,g].begin() + p, p)
					f.begin() + lo1, p,
					g.begin() + lo2, p,
					add_into_iterator<decltype(res.begin()), num>(res.begin() + lo1 + lo2)
				);
				if (i == 2*p-1) break;
				// TODO: Don't recompute if squaring
				multiplier::multiply(
					f.begin() + lo2, p,
					g.begin() + lo1, p,
					add_into_iterator<decltype(res.begin()), num>(res.begin() + lo1 + lo2)
				);
			}
		}
		i++;
	}

	num back() {
		return res[i-1];
	}
};

template <typename num, typename multiplier> struct online_squarer {
	int N; int i;
	std::vector<num> f;
	std::vector<num> res;

	// Computes the first 2N terms of the product
	online_squarer(int N_) : N(N_), i(0), f(N), res(2*N+1, num(0)) {}

	num peek() {
		return res[i];
	}

	void push(num v_f) {
		assert(i < N);
		f[i] = v_f;
		if (i == 0) {
			res[i] += v_f * v_f;
		} else {
			res[i] += 2 * v_f * f[0];
			// TODO: We could do this second half more lazily, since it only affects res[i+1]...
			for (int p = 1; (i & (p-1)) == (p-1); p <<= 1) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				if (i == 2*p-1) {
					multiplier::square(
						// TODO: We can cache FFT([f,g].begin() + p, p)
						f.begin() + lo1, p,
						add_into_iterator<decltype(res.begin()), num>(res.begin() + lo1 + lo2)
					);
					break;
				} else {
					multiplier::multiply(
						// TODO: Use cached FFT
						f.begin() + lo1, p,
						f.begin() + lo2, p,
						add_double_into_iterator<decltype(res.begin()), num>(res.begin() + lo1 + lo2)
					);
				}
			}
		}
		i++;
	}

	num back() {
		return res[i-1];
	}
};

/* namespace ecnerwala */ }
