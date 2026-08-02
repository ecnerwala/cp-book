#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include "fft/engine.hpp"
#include "fft/multiply.hpp"

namespace ecnerwala {

// A polynomial represented by its values evaluated at an Arithmetic Progression (AP).
// TODO: The AP is always assumed to be 0..length-1; store an explicit offset/gap instead?
// Maybe not, this is just more convenient.
template <fft::engine E>
struct ap_sampled_poly : public std::vector<typename E::value_type> {
	using T = typename E::value_type;
	using std::vector<T>::vector;

	int len() const {
		return int(this->size());
	}
	int degree() const {
		return len() - 1;
	}

	ap_sampled_poly& operator += (const ap_sampled_poly& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] += o[i];
		}
		return *this;
	}
	friend ap_sampled_poly operator + (const ap_sampled_poly& a, const ap_sampled_poly& b) {
		assert(a.size() == b.size());
		ap_sampled_poly r(a.size());
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] + b[i];
		}
		return r;
	}
	ap_sampled_poly& operator -= (const ap_sampled_poly& o) {
		assert(len() == o.len());
		for (int i = 0; i < int(o.size()); i++) {
			(*this)[i] -= o[i];
		}
		return *this;
	}
	friend ap_sampled_poly operator - (const ap_sampled_poly& a, const ap_sampled_poly& b) {
		assert(a.size() == b.size());
		ap_sampled_poly r(a.size());
		for (int i = 0; i < r.len(); i++) {
			r[i] = a[i] - b[i];
		}
		return r;
	}

	T eval_at(T k) {
		if (0 <= int(k) && int(k) < len()) {
			return (*this)[int(k)];
		} else {
			// Just do the lagrange interpolation
			std::vector<T> terms(*this);
			{
				// Inverse factorial terms
				T v = 1;
				for (int i = 1; i <= len(); i++) v *= T(i);
				v = inv(v);
				for (int i = len()-1; i >= 0; i--) {
					v *= T(i+1);
					terms[i] *= v;
					terms[len()-1-i] *= (i & 1) ? -v : v;
				}
			}
			{
				// Prefix terms
				T v = 1;
				for (int i = 0; i < len(); i++) {
					terms[i] *= v;
					v *= T(k - i);
				}
			}
			{
				// Suffix terms
				T v = 1;
				for (int i = len() - 1; i >= 0; i--) {
					terms[i] *= v;
					v *= T(k - i);
				}
			}
			T res = 0;
			for (int i = 0; i < len(); i++) res += terms[i];
			return res;
		}
	}

	ap_sampled_poly eval_range(T k, int osz) {
		if (osz == 0) {
			return ap_sampled_poly(osz);
		}
		if (len() == 0) {
			return ap_sampled_poly(osz, T(0));
		}

		// Check for overlaps. We're checking in linear time to avoid unpacking T, but it should be plenty fast.
		// If the field is very very small and we wrap around several times, our runtime can be bad...
		// but then something has already gone wrong, why are you evaluating so many points???
		for (int i = -(len() - 1); i <= osz - 1; i++) {
			if (k+i == T(0)) {
				// everything from [i,i+len()-1) of the output is a match
				ap_sampled_poly res; res.reserve(osz);
				int lo = std::max(0, i);
				int hi = std::min(i+len(), osz);
				{
					auto pref = eval_range(k, lo);
					res.insert(res.end(), pref.begin(), pref.end());
				}
				res.insert(res.end(), this->begin() + (lo - i), this->begin() + (hi - i));
				{
					auto suff = eval_range(k + hi, osz - hi);
					res.insert(res.end(), suff.begin(), suff.end());
				}
				return res;
			}
		}

		std::vector<T> inps(*this);
		{
			// Inverse factorial terms
			T v = 1;
			for (int i = 1; i <= len(); i++) v *= T(i);
			v = inv(v);
			for (int i = len()-1; i >= 0; i--) {
				v *= T(i+1);
				inps[i] *= v;
				inps[len()-1-i] *= (i & 1) ? -v : v;
			}
		}
		std::vector<T> inv_offsets(len() + osz - 1);
		ap_sampled_poly results(osz);
		{
			T v = 1;
			for (int i = - (len() - 1); i <= osz - 1; i++) {
				inv_offsets[i + (len() - 1)] = v;
				v *= k + i;
				if (i >= 0) results[i] = v;
			}
			// Assert there's no overlap
			assert(v != T(0));
			v = inv(v);
			for (int i = osz - 1; i >= -(len() - 1); i--) {
				inv_offsets[i + (len() - 1)] *= v;
				v *= k + i;
				if (i + (len() - 1) <= osz - 1) {
					results[i + (len() - 1)] *= v;
				}
			}
		}
		std::vector<T> prod = fft::middle_product<E>(inv_offsets, inps);
		assert(int(prod.size()) == osz);
		for (int i = 0; i < osz; i++) results[i] *= prod[i];
		return results;
	}

	void extend_right() {
		this->push_back(eval_at(T(len())));
	}
	void extend_left() {
		this->insert(this->begin(), eval_at(T(-1)));
	}

	[[nodiscard]] ap_sampled_poly prefix_sum_inclusive() const {
		ap_sampled_poly r = *this;
		r.extend_right();
		for (int i = 1; i < r.len(); i++) r[i] += r[i-1];
		return r;
	}
};

/* namespace ecnerwala */ }
