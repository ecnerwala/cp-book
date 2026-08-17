#pragma once

#include <algorithm>
#include <cassert>
#include <span>
#include <vector>

#include "fft/multiply.hpp"

namespace wala {

// ==== online multiplication ====

// Online (relaxed) multiplication: computes the first N terms of f*g given the terms one at a time.
template <fft::engine E> struct online_multiplier {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f, g;
	std::vector<T> res;
	std::vector<fft::transformed<E>> f_blocks, g_blocks; // level k: block [2^k, 2^{k+1})

	online_multiplier(int N_) : N(N_), i(0), f(N, T{}), g(N, T{}), res(2*N+1, T{}) {}

	T peek() {
		return res[i];
	}

	void push(T v_f, T v_g) {
		assert(i < N);
		f[i] = v_f;
		g[i] = v_g;
		if (i == 0) {
			res[0] += v_f * v_g;
		} else {
			res[i] += v_f * g[0];
			res[i] += f[0] * v_g;
			for (int p = 1, k = 0; (i & (p-1)) == (p-1); p <<= 1, k++) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				int s = 2*p - 1;
				auto fb = std::span<const T>(f).subspan(p, p);
				auto gb = std::span<const T>(g).subspan(p, p);
				auto out = std::span<T>(res).subspan(lo1 + lo2, s);
				if (i == 2*p-1) {
					f_blocks.emplace_back();
					g_blocks.emplace_back();
					fft::multiply<E>(fb, f_blocks[k], gb, g_blocks[k], out, fft::add_op{});
					break;
				}
				// both products keep f on the left: f_hi * g_lo + f_lo * g_hi
				fft::transformed<E> cf, cg;
				fft::multiply_add2<E>(
						fb, f_blocks[k], std::span<const T>(g).subspan(lo2, p), cg,
						std::span<const T>(f).subspan(lo2, p), cf, gb, g_blocks[k],
						out, fft::add_op{});
			}
		}
		i++;
	}

	T back() {
		return res[i-1];
	}
};

template <fft::engine E> struct online_squarer {
	using T = typename E::value_type;
	int N; int i;
	std::vector<T> f;
	std::vector<T> res;
	std::vector<fft::transformed<E>> f_blocks;

	online_squarer(int N_) : N(N_), i(0), f(N, T{}), res(2*N+1, T{}) {}

	T peek() {
		return res[i];
	}

	void push(T v_f) {
		assert(i < N);
		f[i] = v_f;
		if (i == 0) {
			res[0] += v_f * v_f;
		} else {
			if constexpr (E::commutative) res[i] += (v_f + v_f) * f[0];
			else res[i] += v_f * f[0] + f[0] * v_f;
			for (int p = 1, k = 0; (i & (p-1)) == (p-1); p <<= 1, k++) {
				int lo1 = p;
				int lo2 = i + 1 - p;
				int s = 2*p - 1;
				auto fb = std::span<const T>(f).subspan(p, p);
				auto fw = std::span<const T>(f).subspan(lo2, p);
				auto out = std::span<T>(res).subspan(lo1 + lo2, s);
				if (i == 2*p-1) {
					f_blocks.emplace_back();
					fft::square<E>(fb, f_blocks[k], out, fft::add_op{});
					break;
				}
				fft::transformed<E> cw;
				if constexpr (E::commutative) {
					fft::multiply<E>(fb, f_blocks[k], fw, cw, out, fft::add_twice_op{});
				} else {
					// f_hi * f_lo + f_lo * f_hi from the same two transforms
					fft::multiply_add2<E>(fb, f_blocks[k], fw, cw,
							fw, cw, fb, f_blocks[k], out, fft::add_op{});
				}
			}
		}
		i++;
	}

	T back() {
		return res[i-1];
	}
};

/* namespace wala */ }
