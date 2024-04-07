#pragma once

/*
 * This is mostly inspired by https://golang.org/src/index/suffixarray/sais.go.
 */

#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include <type_traits>

#include "rmq.hpp"

template<class T> int sz(T&& arg) { using std::size; return int(size(std::forward<T>(arg))); }

class SuffixArray {
public:
	using index_t = int;
	int N;
	std::vector<index_t> sa;
	std::vector<index_t> rank;
	// lcp[i] = get_lcp(sa[i], sa[i+1])
	std::vector<index_t> lcp;
	RangeMinQuery<std::pair<index_t, index_t>> rmq;

	SuffixArray() {}

	template <typename String> static SuffixArray construct_raw(const String& S, index_t sigma) {
		int N = sz(S);
		SuffixArray sa(N);

		sa.build_sa(S, sigma);
		sa.build_rank();
		sa.build_lcp(S);
		sa.build_rmq();

		return sa;
	}

	// Pass a function which returns a value in [0, sigma)
	template <typename String, typename F> static SuffixArray map_and_construct(const String& S, const F& f, int sigma) {
		std::vector<decltype((f(S[0])))> mapped(sz(S));
		for (int i = 0; i < sz(S); i++) {
			mapped[i] = f(S[i]);
			assert(0 <= int(mapped[i]) && int(mapped[i]) < sigma);
		}
		return construct_raw(mapped, sigma);
	}

	// Sorts the elements of S and then runs suffix array. This takes O(N log N) time with no dependence on sigma.
	template <typename String> static SuffixArray sort_and_construct(const String& S) {
		using std::begin;
		using std::end;
		using value_type = typename std::iterator_traits<decltype(begin(S))>::value_type;
		using compressed_value_type = typename std::conditional<
			sizeof(value_type) < sizeof(index_t),
			value_type,
			index_t
		>::type;

		std::vector<compressed_value_type> compressed_s(sz(S));
		int sigma = 0;

		{
			std::vector<value_type> vals(begin(S), end(S));
			std::sort(vals.begin(), vals.end());
			vals.resize(unique(vals.begin(), vals.end()) - vals.begin());
			for (int i = 0; i < sz(S); i++) {
				compressed_s[i] = compressed_value_type(index_t(std::lower_bound(vals.begin(), vals.end(), S[i]) - vals.begin()));
			}
			sigma = int(vals.size());
		}

		return construct_raw(compressed_s, sigma);
	}

	// Shifts the elements so that sigma = max(S) - min(S) + 1
	template <typename String> static SuffixArray shift_and_construct(const String& S) {
		using std::begin;
		using std::end;
		using value_type = typename std::iterator_traits<decltype(begin(S))>::value_type;

		std::vector<value_type> compressed_s(sz(S));
		int sigma = 0;

		if (sz(S) > 0) {
			value_type lo = *begin(S), hi = *begin(S);
			for (const auto& x : S) {
				if (x < lo) lo = x;
				if (x > hi) hi = x;
			}

			for (int i = 0; i < sz(S); i++) {
				compressed_s[i] = value_type(S[i] - lo);
			}
			sigma = int(hi - lo + 1);
		}

		return construct_raw(compressed_s, sigma);
	}

	// Renumber/filter to only the used elements with bucket sorting. Still takes O(max(S) - min(S) + 1) memory/time,
	// but should be less memory than `shift_and_construct` when sigma ~ N and max(S) - min(S) + 1 > N.
	template <typename String> static SuffixArray bucket_and_construct(const String& S) {
		using std::begin;
		using std::end;
		using value_type = typename std::iterator_traits<decltype(begin(S))>::value_type;
		using compressed_value_type = typename std::conditional<
			sizeof(value_type) < sizeof(index_t),
			value_type,
			index_t
		>::type;

		std::vector<compressed_value_type> compressed_s(sz(S));
		int sigma = 0;

		if (sz(S) > 0) {
			value_type lo = *begin(S), hi = *begin(S);
			for (const auto& x : S) {
				if (x < lo) lo = x;
				if (x > hi) hi = x;
			}

			std::vector<compressed_value_type> buckets(hi - lo + 1, 0);
			for (const auto& x : S) {
				buckets[x - lo] = 1;
			}
			for (int v = 0; v < int(buckets.size()); v++) {
				if (buckets[v]) buckets[v] = compressed_value_type(sigma++);
			}

			for (int i = 0; i < sz(S); i++) {
				compressed_s[i] = buckets[S[i] - lo];
			}
		}

		return construct_raw(compressed_s, sigma);
	}

	index_t get_lcp(index_t a, index_t b) const {
		if (a == b) return N-a;
		a = rank[a], b = rank[b];
		if (a > b) std::swap(a, b);
		return rmq.query(a, b-1).first;
	}

	// Get the split in the suffix tree, using half-open intervals
	// Returns len, idx
	std::pair<index_t, index_t> get_split(index_t l, index_t r) const {
		assert(r - l > 1);
		return rmq.query(l, r-2);
	}

private:
	explicit SuffixArray(int N_) : N(N_) {}

	template <typename String> void build_sa(const String& S, index_t sigma) {
		sa = std::vector<index_t>(N+1);
		assert(sigma >= 0);
		for (auto s : S) assert(0 <= index_t(s) && index_t(s) < sigma);
		std::vector<index_t> tmp(sigma + std::max(N, sigma));
		SuffixArray::sais<String>(N, S, sa.data(), sigma, tmp.data());
	}

	template <typename String> static void sais(int N, const String& S, index_t* sa, int sigma, index_t* tmp) {
		if (N == 0) {
			sa[0] = 0;
			return;
		} else if (N == 1) {
			sa[0] = 1;
			sa[1] = 0;
			return;
		}

		// Phase 1: Initialize the frequency array, which will let us lookup buckets.
		index_t* freq = tmp; tmp += sigma;
		memset(freq, 0, sizeof(*freq) * sigma);
		for (int i = 0; i < N; i++) {
			++freq[index_t(S[i])];
		}
		auto build_bucket_start = [&]() {
			int cur = 1;
			for (int v = 0; v < sigma; v++) {
				tmp[v] = cur;
				cur += freq[v];
			}
		};
		auto build_bucket_end = [&]() {
			int cur = 1;
			for (int v = 0; v < sigma; v++) {
				cur += freq[v];
				tmp[v] = cur;
			}
		};

		int num_pieces = 0;

		int first_endpoint = 0;
		// Phase 2: find the right-endpoints of the pieces
		{
			build_bucket_end();

			// Initialize the final endpoint out-of-band this way so that we don't try to look up tmp[-1].
			// This doesn't count towards num_pieces.
			sa[0] = N;

			index_t c0 = S[N-1], c1 = -1; bool isS = false;
			for (int i = N-2; i >= 0; i--) {
				c1 = c0;
				c0 = S[i];
				if (c0 < c1) {
					isS = true;
				} else if (c0 > c1 && isS) {
					isS = false;
					// insert i+1
					sa[first_endpoint = --tmp[c1]] = i+1;
					++num_pieces;
				}
			}
		}

		// If num_pieces <= 1, we don't need to actually run the recursion, it's just sorted automatically
		// Otherwise, we're going to rebucket
		if (num_pieces > 1) {
			// Remove the first endpoint, we don't need to run the IS on this
			sa[first_endpoint] = 0;

			// Run IS for L-type
			{
				build_bucket_start();
				for (int z = 0; z <= N; z++) {
					int v = sa[z];
					if (!v) continue;

					// Leave for the S-round
					if (v < 0) continue;

					// clear out our garbage
					sa[z] = 0;

					--v;
					index_t c0 = S[v-1], c1 = S[v];
					sa[tmp[c1]++] = (c0 < c1) ? ~v : v;
				}
			}

			index_t* const sa_end = sa + N + 1;

			index_t* pieces = sa_end;
			// Run IS for S-type and compactify
			{
				build_bucket_end();
				for (int z = N; z >= 0; z--) {
					int v = sa[z];
					if (!v) continue;

					// clear our garbage
					sa[z] = 0;

					if (v > 0) {
						*--pieces = v;
						continue;
					}

					v = ~v;

					--v;
					index_t c0 = S[v-1], c1 = S[v];
					sa[--tmp[c1]] = (c0 > c1) ? v : ~v;
				}
			}

			// Compute the lengths of the pieces in preparation for equality
			// comparison, and store them in sa[v/2]. We set the length of the
			// final piece to 0; it compares unequal to everything because of
			// the sentinel.
			{
				int prv_start = N;
				index_t c0 = S[N-1], c1 = -1; bool isS = false;
				for (int i = N-2; i >= 0; i--) {
					c1 = c0;
					c0 = S[i];
					if (c0 < c1) {
						isS = true;
					} else if (c0 > c1 && isS) {
						isS = false;

						// insert i+1
						int v = i+1;
						sa[v>>1] = prv_start == N ? 0 : prv_start - v;
						prv_start = v;
					}
				}
			}

			// Compute the alphabet, storing the result into sa[v/2].
			int next_sigma = 0;
			{
				int prv_len = -1, prv_v = 0;
				for (int i = 0; i < num_pieces; i++) {
					int v = pieces[i];
					int len = sa[v>>1];

					bool eq = prv_len == len;
					for (int a = 0; eq && a < len; ++a) {
						eq = S[v+a] == S[prv_v+a];
					}
					if (!eq) {
						next_sigma++;
						prv_len = len;
						prv_v = v;
					}

					sa[v>>1] = next_sigma; // purposely leave this 1 large to check != 0
				}
			}

			if (next_sigma == num_pieces) {
				sa[0] = N;
				memcpy(sa+1, pieces, sizeof(*sa) * num_pieces);
			} else {
				index_t* next_S = sa_end;

				// Finally, pack the input to the SA
				{
					for (int i = (N-1)>>1; i >= 0; i--) {
						int v = sa[i];
						if (v) *--next_S = v-1;
						sa[i] = 0;
					}
				}

				memset(sa, 0, sizeof(*sa) * (num_pieces+1));
				sais<const index_t*>(num_pieces, next_S, sa, next_sigma, tmp);

				{ // Compute the piece start points again and use those to map up the suffix array
					next_S = sa_end;
					index_t c0 = S[N-1], c1 = -1; bool isS = false;
					for (int i = N-2; i >= 0; i--) {
						c1 = c0;
						c0 = S[i];
						if (c0 < c1) {
							isS = true;
						} else if (c0 > c1 && isS) {
							isS = false;

							int v = i+1;
							*--next_S = v;
						}
					}
					sa[0] = N;
					for (int i = 1; i <= num_pieces; i++) {
						sa[i] = next_S[sa[i]];
					}
				}
			}

			// zero everything else
			memset(sa+num_pieces+1, 0, sizeof(*sa) * (N - num_pieces));

			{
				// Scatter the finished pieces
				build_bucket_end();
				for (int i = num_pieces; i > 0; i--) {
					int v = sa[i];
					sa[i] = 0;

					index_t c1 = S[v];
					sa[--tmp[c1]] = v;
				}
			}
		}

		// Home stretch! Just finish out with the L-type and then S-type
		{
			build_bucket_start();
			for (int z = 0; z <= N; z++) {
				int v = sa[z];
				if (v <= 0) continue;
				--v;
				index_t c1 = S[v];
				index_t c0 = v ? S[v-1] : c1; // if v = 0, we don't want to invert
				sa[tmp[c1]++] = (c0 < c1) ? ~v : v;
			}
		}

		// This just aggressively overwrites our original scattered pieces with the correct values
		{
			build_bucket_end();
			for (int z = N; z >= 0; z--) {
				int v = sa[z];
				if (v >= 0) continue;
				sa[z] = v = ~v;
				--v;
				index_t c1 = S[v];
				index_t c0 = v ? S[v-1] : c1+1;
				sa[--tmp[c1]] = (c0 > c1) ? v : ~v;
			}
		}
	}

	void build_rank() {
		rank = std::vector<index_t>(N+1);
		for (int i = 0; i <= N; i++) rank[sa[i]] = i;
	}

	template <typename String> void build_lcp(const String& S) {
		assert(sz(S) == N);
		lcp = std::vector<index_t>(N);
		for (int i = 0, k = 0; i < N - 1; i++) {
			int j = sa[rank[i]-1];
			while (k < N - std::max(i, j) && S[i+k] == S[j+k]) k++;
			lcp[rank[i]-1] = k;
			if (k) --k;
		}
	}

	void build_rmq() {
		std::vector<std::pair<index_t, index_t>> lcp_idx(N);
		for (int i = 0; i < N; i++) {
			lcp_idx[i] = {lcp[i], i+1};
		}
		rmq = RangeMinQuery<std::pair<index_t, index_t>>(std::move(lcp_idx));
	}
};

class PrefixArray : private SuffixArray {
	PrefixArray(const SuffixArray& sa_) : SuffixArray(sa_) {}
	PrefixArray(SuffixArray&& sa_) : SuffixArray(std::move(sa_)) {}
public:
	PrefixArray() {}
	template <typename String> static PrefixArray construct_raw(const String& S, int sigma) {
		return PrefixArray(SuffixArray::construct_raw(String(S.rbegin(), S.rend())), sigma);
	}

	// TODO: Fill in other constructors

	int get_lcs(int a, int b) const {
		return SuffixArray::get_lcp(SuffixArray::N - a, SuffixArray::N - b);
	}
};
