#pragma once

#include <algorithm>
#include <iterator>
#include <span>
#include <utility>
#include <vector>

/**
 * Author: Andrew He
 * Source: http://neerc.ifmo.ru/trains/toulouse/2017/fft2.pdf
 * Papers about accuracy: http://www.daemonology.net/papers/fft.pdf, http://www.cs.berkeley.edu/~fateman/papers/fftvsothers.pdf
 * For integers rounding works if $(|a| + |b|)\max(a, b) < \mathtt{\sim} 10^9$, or in theory maybe $10^6$.
 *
 * Abstraction layers:
 *   fft_core<num>   FFT itself and other ops on rings with 2^k-th roots of unity. We use bit-reversed indexing in the frequency domain.
 *
 *   engines         Engines for packing/unpacking arbitrary rings for convolution (the `engine` concept).
 *                   Still expose (opaque) transform-domain objects for caching/fusion.
 *
 *   multiply layer  Wrappers for convolving bounded sequences: track length/truncation.
 *
 *   value types     series::vec<E, exact> - R[[x]]
 *                   series::exact<E> - exact (finite-support) power series
 *                   series::trunc<E> - truncated prefix of an (infinite) power series
 *
 *                   polynomials - R[x]. Under x -> 1/x a polynomial becomes a Laurent polynomial in 1/x;
 *                                 shifting by x^{deg P} (reversal) lands it in R[[x]], and we store that exact series.
 *                   poly::vec<E> - polynomial type, supporting natural indexing
 *                   poly::form<E> - finite-support linear forms, via the pairing <P, S> = [x^0] P(1/x) S(x)
 *                                    a linear form is one side of this pairing, applied to the other
 *
 *                   online_multiplier<E> - online (relaxed) multiplication of 2 sequences in n log^2 n time
 *                   ap_sampled_poly<E> - a polynomial stored as its evaluations on an arithmetic progression
 */

namespace ecnerwala {

template <class T> int sz(T&& arg) {
	using std::size;
	return int(size(std::forward<T>(arg)));
}
inline int nextPow2(int s) { return 1 << (s > 1 ? 32 - __builtin_clz(s - 1) : 0); }

namespace fft {

using std::swap;
using std::vector;
using std::min;
using std::max;

// Reusable scratch buffers. Not thread-safe by default: this is deliberately plain
// static storage so single-threaded programs pay no TLS indirection; define
// ECNERWALA_FFT_POOL_STORAGE to `thread_local` for multithreaded use.
#ifndef ECNERWALA_FFT_POOL_STORAGE
#define ECNERWALA_FFT_POOL_STORAGE
#endif
template <typename T> struct buffer_pool {
	static inline ECNERWALA_FFT_POOL_STORAGE std::vector<std::vector<T>> free_list;
	struct handle {
		std::vector<T> v;
		explicit handle(int n) {
			if (!free_list.empty()) {
				v = std::move(free_list.back());
				free_list.pop_back();
			}
			v.assign(n, T());
		}
		handle(const handle&) = delete;
		handle& operator = (const handle&) = delete;
		handle(handle&& o) noexcept : v(std::move(o.v)) {}
		~handle() {
			if (v.capacity()) free_list.push_back(std::move(v));
		}
		T& operator [] (int i) { return v[i]; }
		operator std::span<T>() { return std::span<T>(v); }
		std::span<T> span() { return std::span<T>(v); }
	};
	static handle get(int n) { return handle(n); }
};

/* namespace fft */ }

/* namespace ecnerwala */ }
