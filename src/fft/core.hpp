#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <span>
#include <vector>
#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

#include "fft/common.hpp"
#include "modnum.hpp"

namespace ecnerwala::fft {

// ==== core: roots, buffers, raw transforms ====

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
// 998244353 = (119 << 23) + 1 = 2^30 - 2^26 - 2^23 + 1
template <> struct primitive_root<998244353> {
	static const int value = 3;
};
// babybear prime
template <> struct primitive_root<(15 << 27) + 1> {
	static const int value = 31;
};
// koalabear prime
template <> struct primitive_root<(127 << 24) + 1> {
	static const int value = 3;
};
template <> struct primitive_root<(7 << 26) + 1> {
	static const int value = 3;
};
template <> struct primitive_root<(5 << 25) + 1> {
	static const int value = 3;
};
template <int MOD> struct getRoot<modnum<MOD>> {
	static modnum<MOD> f(int k) {
		assert((MOD-1)%k == 0);
		return power(modnum<MOD>(primitive_root<MOD>::value), (MOD-1)/k);
	}
};
template <> struct getRoot<mod_goldilocks> {
	static mod_goldilocks f(int k) {
		assert((mod_goldilocks::MOD-1)%k == 0);
		return power(mod_goldilocks(mod_goldilocks::PRIMITIVE_ROOT), (mod_goldilocks::MOD-1)/k);
	}
};

// SIMD (AVX2) butterflies for 32-bit prime moduli: twiddle multiplies use Shoup's
// precomputed-quotient trick (w' = floor(w 2^32 / MOD)), so the data stays in the
// normal representation, lazily in [0, 2 MOD) between passes and normalized back
// to [0, MOD) on the final pass. Requires MOD < 2^30 so lazy values fit 32 bits.
// The SIMD functions carry target("avx2") and are gated at runtime, so no special
// compile flags are needed.
template <typename num> constexpr bool ntt_simd = false;
#if defined(__x86_64__) || defined(__i386__)
template <int MOD> constexpr bool ntt_simd<modnum<MOD>> = MOD < (1 << 30);
inline const bool cpu_has_avx2 = __builtin_cpu_supports("avx2");
#else
inline const bool cpu_has_avx2 = false;
#endif

// We take the bit-reverse convention: the coefficient of a[i] -> b[j] is omega^{i * bit_reverse(j)}.
// This means that the size 2^{k-1} transform is the prefix of the size 2^k transform (wrapping the input).
//
// We mostly work with spans here:
//   spans of transforms are expected to have length exactly 2^k
//   spans of inputs/outputs are expected to have length [0, 2^{k+1})
//
// Inputs/outputs are treated mod x^{2^k} - 1.
// Their length is allowed to be bigger than 2^k mostly to perform ops on sequences of size n+1 with only transforms of size n.
// The upper bound of 2^{k+1} is arbitrary: we could tighten it to 2^k + 1 or loosen it to infinity, this is just a "defensive" choice.
template <typename num> struct fft_core {

	static inline vector<int> rev;
	// rt[2^k + i] = 1^{i / 2^(k+1)}
	// TODO: can we get rid of inv_rt; alternatively, should we store inv_rt in bit-reverse order?
	static inline vector<num> rt, inv_rt;
	static inline vector<uint32_t> rt_sh, inv_rt_sh;

	static void init(int n) {
		if (n <= sz(rt)) return;
		rev.resize(n);
		for (int i = 0; i < n; i++) {
			rev[i] = (rev[i>>1] | ((i&1)*n)) >> 1;
		}
		rt.reserve(n); inv_rt.reserve(n);
		while (sz(rt) < 2 && sz(rt) < n) { rt.push_back(num(1)); inv_rt.push_back(num(1)); }
		for (int k = sz(rt); k < n; k *= 2) {
			rt.resize(2*k); inv_rt.resize(2*k);
			num z = getRoot<num>::f(2*k);
			num iz = inv(z);
			for (int i = k/2; i < k; i++) {
				rt[2*i] = rt[i], rt[2*i+1] = rt[i]*z;
				inv_rt[2*i] = inv_rt[i], inv_rt[2*i+1] = inv_rt[i]*iz;
			}
		}
		if constexpr (ntt_simd<num>) {
			rt_sh.resize(size_t(sz(rt))); inv_rt_sh.resize(size_t(sz(inv_rt)));
			for (int i = 0; i < sz(rt); i++) {
				rt_sh[size_t(i)] = uint32_t((uint64_t(uint32_t(int(rt[i]))) << 32) / uint32_t(num::MOD));
				inv_rt_sh[size_t(i)] = uint32_t((uint64_t(uint32_t(int(inv_rt[i]))) << 32) / uint32_t(num::MOD));
			}
		}
	}

	// bit-reversal of i as a log2(n)-bit number
	static int brev(int i, int n) {
		int s = __builtin_ctz(unsigned(sz(rev)/n));
		return rev[i] >> s;
	}
	// index of the conjugate evaluation point, in the returned bit-reversed order
	static int conj_index(int j) {
		return j == 0 ? 0 : j ^ ((1 << (31 - __builtin_clz(unsigned(j)))) - 1);
	}

#if defined(__x86_64__) || defined(__i386__)
	// q = mulhi32(a, b) per lane
	[[gnu::target("avx2")]] static __m256i mulhi_u32(__m256i a, __m256i b) {
		__m256i lo = _mm256_srli_epi64(_mm256_mul_epu32(a, b), 32);
		__m256i hi = _mm256_mul_epu32(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32));
		return _mm256_blend_epi32(lo, hi, 0b10101010);
	}
	// a * w mod MOD via Shoup, result in [0, 2 MOD); requires a < 2^32, w < MOD
	[[gnu::target("avx2")]] static __m256i mul_shoup(__m256i a, __m256i w, __m256i wsh, __m256i vP) {
		__m256i q = mulhi_u32(a, wsh);
		return _mm256_sub_epi32(_mm256_mullo_epi32(a, w), _mm256_mullo_epi32(q, vP));
	}
	// [0, 2 lim) -> [0, lim)
	[[gnu::target("avx2")]] static __m256i lazy_reduce(__m256i x, __m256i vlim) {
		return _mm256_min_epu32(x, _mm256_sub_epi32(x, vlim));
	}
	// out[i] = a[i] * w mod MOD for a constant w, inputs in [0, MOD), outputs in [0, MOD)
	[[gnu::target("avx2")]] static void scale_simd(const num* a_, num* out_, int n, num w) {
		const uint32_t* a = reinterpret_cast<const uint32_t*>(a_);
		uint32_t* out = reinterpret_cast<uint32_t*>(out_);
		const uint32_t P = uint32_t(num::MOD);
		uint32_t wv = uint32_t(int(w));
		uint32_t wsh = uint32_t((uint64_t(wv) << 32) / P);
		const __m256i vP = _mm256_set1_epi32(int(P));
		const __m256i vw = _mm256_set1_epi32(int(wv)), vwsh = _mm256_set1_epi32(int(wsh));
		int i = 0;
		for (; i + 8 <= n; i += 8) {
			__m256i x = _mm256_loadu_si256((const __m256i*)(a + i));
			__m256i r = lazy_reduce(mul_shoup(x, vw, vwsh, vP), vP);
			_mm256_storeu_si256((__m256i*)(out + i), r);
		}
		for (; i < n; i++) {
			uint32_t q = uint32_t((uint64_t(a[i]) * wsh) >> 32);
			uint32_t r = a[i] * wv - q * P;
			out[i] = r >= P ? r - P : r;
		}
	}
	// b[i] = coeffs[i] * rt[n + i], inputs in [0, MOD), outputs in [0, MOD)
	[[gnu::target("avx2")]] static void twiddle_simd(const num* c_, num* b_, int n0, int lo) {
		const uint32_t* c = reinterpret_cast<const uint32_t*>(c_);
		uint32_t* b = reinterpret_cast<uint32_t*>(b_);
		const uint32_t* w_ = reinterpret_cast<const uint32_t*>(rt.data());
		const uint32_t P = uint32_t(num::MOD);
		const __m256i vP = _mm256_set1_epi32(int(P));
		int i = 0;
		for (; i + 8 <= lo; i += 8) {
			__m256i x = _mm256_loadu_si256((const __m256i*)(c + i));
			__m256i w = _mm256_loadu_si256((const __m256i*)(w_ + n0 + i));
			__m256i wsh = _mm256_loadu_si256((const __m256i*)(rt_sh.data() + n0 + i));
			__m256i r = lazy_reduce(mul_shoup(x, w, wsh, vP), vP);
			_mm256_storeu_si256((__m256i*)(b + i), r);
		}
		for (; i < lo; i++) {
			uint32_t q = uint32_t((uint64_t(c[i]) * rt_sh[size_t(n0 + i)]) >> 32);
			uint32_t r = c[i] * w_[n0 + i] - q * P;
			b[i] = r >= P ? r - P : r;
		}
	}
	[[gnu::target("avx2")]] static void forward_simd(std::span<num> a_) {
		int n = sz(a_);
		uint32_t* a = reinterpret_cast<uint32_t*>(a_.data());
		const uint32_t P = uint32_t(num::MOD);
		const __m256i vP = _mm256_set1_epi32(int(P)), v2P = _mm256_set1_epi32(int(2 * P));
		const uint32_t* w_ = reinterpret_cast<const uint32_t*>(rt.data());
		for (int k = n/2; k >= 8; k /= 2) {
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j += 8) {
					__m256i u = _mm256_loadu_si256((const __m256i*)(a + i + j));
					__m256i v = _mm256_loadu_si256((const __m256i*)(a + i + j + k));
					__m256i w = _mm256_loadu_si256((const __m256i*)(w_ + j + k));
					__m256i wsh = _mm256_loadu_si256((const __m256i*)(rt_sh.data() + j + k));
					__m256i s = lazy_reduce(_mm256_add_epi32(u, v), v2P);
					__m256i d = _mm256_add_epi32(_mm256_sub_epi32(u, v), v2P);
					_mm256_storeu_si256((__m256i*)(a + i + j), s);
					_mm256_storeu_si256((__m256i*)(a + i + j + k), mul_shoup(d, w, wsh, vP));
				}
			}
		}
		for (int k = min(n/2, 4); k >= 1; k /= 2) {
			bool last = (k == 1);
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j++) {
					uint32_t u = a[i+j], v = a[i+j+k];
					uint32_t s = u + v; if (s >= 2*P) s -= 2*P;
					uint32_t d = u - v + 2*P;
					uint32_t q = uint32_t((uint64_t(d) * rt_sh[size_t(j+k)]) >> 32);
					uint32_t t = d * w_[j+k] - q * P;
					if (last) { if (s >= P) s -= P; if (t >= P) t -= P; }
					a[i+j] = s;
					a[i+j+k] = t;
				}
			}
		}
	}
	[[gnu::target("avx2")]] static void inverse_simd(std::span<num> a_) {
		int n = sz(a_);
		uint32_t* a = reinterpret_cast<uint32_t*>(a_.data());
		const uint32_t P = uint32_t(num::MOD);
		const __m256i vP = _mm256_set1_epi32(int(P)), v2P = _mm256_set1_epi32(int(2 * P));
		const uint32_t* w_ = reinterpret_cast<const uint32_t*>(inv_rt.data());
		for (int k = 1; k < n && k < 8; k *= 2) {
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j++) {
					uint32_t y = a[i+j+k];
					uint32_t q = uint32_t((uint64_t(y) * inv_rt_sh[size_t(j+k)]) >> 32);
					uint32_t t = y * w_[j+k] - q * P;
					uint32_t x = a[i+j];
					uint32_t s = x + t; if (s >= 2*P) s -= 2*P;
					uint32_t d = x - t + 2*P; if (d >= 2*P) d -= 2*P;
					a[i+j] = s;
					a[i+j+k] = d;
				}
			}
		}
		for (int k = 8; k < n; k *= 2) {
			bool last = (2*k == n);
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j += 8) {
					__m256i x = _mm256_loadu_si256((const __m256i*)(a + i + j));
					__m256i y = _mm256_loadu_si256((const __m256i*)(a + i + j + k));
					__m256i w = _mm256_loadu_si256((const __m256i*)(w_ + j + k));
					__m256i wsh = _mm256_loadu_si256((const __m256i*)(inv_rt_sh.data() + j + k));
					__m256i t = mul_shoup(y, w, wsh, vP);
					__m256i s = lazy_reduce(_mm256_add_epi32(x, t), v2P);
					__m256i d = lazy_reduce(_mm256_add_epi32(_mm256_sub_epi32(x, t), v2P), v2P);
					if (last) { s = lazy_reduce(s, vP); d = lazy_reduce(d, vP); }
					_mm256_storeu_si256((__m256i*)(a + i + j), s);
					_mm256_storeu_si256((__m256i*)(a + i + j + k), d);
				}
			}
		}
		if (n < 16) {
			for (int i = 0; i < n; i++) { if (a[i] >= P) a[i] -= P; }
		}
	}
#endif

	static void forward(std::span<num> a) {
		int n = sz(a);
		if (n <= 1) return;
		init(n);
		if constexpr (ntt_simd<num>) { if (cpu_has_avx2) { forward_simd(a); return; } }
		for (int k = n/2; k >= 1; k /= 2) {
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j++) {
					num u = a[i+j], v = a[i+j+k];
					a[i+j] = u + v;
					a[i+j+k] = (u - v) * rt[j+k];
				}
			}
		}
	}

	static void inverse(std::span<num> a) {
		int n = sz(a);
		if (n <= 1) return;
		init(n);
		if constexpr (ntt_simd<num>) { if (cpu_has_avx2) { inverse_simd(a); return; } }
		for (int k = 1; k < n; k *= 2) {
			for (int i = 0; i < n; i += 2*k) {
				for (int j = 0; j < k; j++) {
					num t = inv_rt[j+k] * a[i+j+k];
					a[i+j+k] = a[i+j] - t;
					a[i+j] = a[i+j] + t;
				}
			}
		}
	}

	// Extend a size 2^{k-1} transform to size 2^k; we need the coefficients.
	// t must have size 2^k, and coeffs must have size at most 2^{k+1}.
	static void extend(std::span<num> t, std::span<const num> coeffs) {
		int n = sz(t) / 2;
		assert(sz(coeffs) <= 2 * n);
		init(sz(t));
		auto b = t.subspan(n, n);
		int lo = min(sz(coeffs), n);
		bool simd = false;
		if constexpr (ntt_simd<num>) {
			if (cpu_has_avx2) { twiddle_simd(coeffs.data(), b.data(), n, lo); simd = true; }
		}
		if (!simd) {
			for (int i = 0; i < lo; i++) {
				// rt[n + i] = w_{2n}^i
				b[i] = coeffs[i] * rt[n + i];
			}
		}
		std::fill(b.begin() + lo, b.end(), num(0));
		for (int i = n; i < sz(coeffs); i++) {
			b[i - n] = b[i - n] - coeffs[i] * rt[i];
		}
		forward(b);
	}

	// Consider t = transform(P) and P(x) = E(x^2) + O(x^2) * x
	// `even_half` and `odd_half` extract a size 2^{k-1} transform of E/O, respectively.
	static void even_half(std::span<const num> t, std::span<num> out) {
		int n = sz(out);
		assert(sz(t) >= 2*n);
		num half = inv(num(2));
		for (int j = 0; j < n; j++) out[j] = (t[2*j] + t[2*j+1]) * half;
	}
	static void odd_half(std::span<const num> t, std::span<num> out) {
		int n = sz(out);
		assert(sz(t) >= 2*n);
		init(2*n);
		num half = inv(num(2));
		for (int j = 0; j < n; j++) {
			// entry j of the size-2n transform pairs (w, -w) with w = w_{2n}^{brev(j, n)}
			out[j] = (t[2*j] - t[2*j+1]) * half * inv_rt[n + brev(j, n)];
		}
	}
};

/* namespace ecnerwala::fft */ }
