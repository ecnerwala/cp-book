#pragma once

/*
 * This is mostly inspired by https://golang.org/src/index/suffixarray/sais.go.
 */

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstring>
#include <iterator>
#include <ranges>
#include <type_traits>

#include "containers/vec.hpp"

#include "ds/rmq.hpp"

namespace wala {

template <typename S> concept text = std::ranges::random_access_range<S> && std::ranges::sized_range<S>;
template <typename S> concept integer_text = text<S> && std::integral<std::ranges::range_value_t<S>>;

// Layered suffix array: SuffixArrayBase computes just sa/rank, each further
// layer statically opts into one more derived structure. Use the leaf classes
// SuffixArray, SuffixArrayLCP, or SuffixArrayRMQ; the named constructors on
// each return that type.
template <typename Self> class SuffixArrayBase {
public:
	using index_t = int;
	int N;
	vec<index_t> sa;
	vec<index_t> rank;

	SuffixArrayBase() : N(0) {}

	template <integer_text String> static Self construct_raw(const String& S, index_t sigma) {
		Self res;
		res.build(S, sigma);
		return res;
	}

	// Pass a function which returns a value in [0, sigma)
	template <text String, std::regular_invocable<std::ranges::range_reference_t<const String>> F> static Self map_and_construct(const String& S, const F& f, int sigma) {
		auto s = std::ranges::begin(S);
		vec<std::invoke_result_t<const F&, std::ranges::range_reference_t<const String>>> mapped(int(std::ranges::size(S)), uninit);
		for (int i = 0; i < int(std::ranges::size(S)); i++) {
			mapped[i] = f(s[i]);
			assert(0 <= int(mapped[i]) && int(mapped[i]) < sigma);
		}
		return construct_raw(mapped, sigma);
	}

	// Sorts the elements of S and then runs suffix array. This takes O(N log N) time with no dependence on sigma.
	template <text String> static Self sort_and_construct(const String& S) requires std::totally_ordered<std::ranges::range_value_t<String>> {
		auto s = std::ranges::begin(S);
		using value_type = std::ranges::range_value_t<String>;
		using compressed_value_type = std::conditional<
			sizeof(value_type) < sizeof(index_t),
			value_type,
			index_t
		>::type;

		vec<compressed_value_type> compressed_s(int(std::ranges::size(S)), uninit);
		int sigma = 0;

		{
			vec<value_type> vals(std::from_range, S);
			std::sort(vals.begin(), vals.end());
			auto vals_end = std::unique(vals.begin(), vals.end());
			for (int i = 0; i < int(std::ranges::size(S)); i++) {
				compressed_s[i] = compressed_value_type(index_t(std::lower_bound(vals.begin(), vals_end, s[i]) - vals.begin()));
			}
			sigma = int(vals_end - vals.begin());
		}

		return construct_raw(compressed_s, sigma);
	}

	// Shifts the elements so that sigma = max(S) - min(S) + 1
	template <integer_text String> static Self shift_and_construct(const String& S) {
		auto s = std::ranges::begin(S);
		using value_type = std::ranges::range_value_t<String>;

		vec<value_type> compressed_s(int(std::ranges::size(S)), uninit);
		int sigma = 0;

		if (int(std::ranges::size(S)) > 0) {
			value_type lo = s[0], hi = s[0];
			for (const auto& x : S) {
				if (x < lo) lo = x;
				if (x > hi) hi = x;
			}

			for (int i = 0; i < int(std::ranges::size(S)); i++) {
				compressed_s[i] = value_type(s[i] - lo);
			}
			sigma = int(hi - lo + 1);
		}

		return construct_raw(compressed_s, sigma);
	}

	// Renumber/filter to only the used elements with bucket sorting. Still takes O(max(S) - min(S) + 1) memory/time,
	// but should be less memory than `shift_and_construct` when sigma ~ N and max(S) - min(S) + 1 > N.
	template <integer_text String> static Self bucket_and_construct(const String& S) {
		auto s = std::ranges::begin(S);
		using value_type = std::ranges::range_value_t<String>;
		using compressed_value_type = std::conditional<
			sizeof(value_type) < sizeof(index_t),
			value_type,
			index_t
		>::type;

		vec<compressed_value_type> compressed_s(int(std::ranges::size(S)), uninit);
		int sigma = 0;

		if (int(std::ranges::size(S)) > 0) {
			value_type lo = s[0], hi = s[0];
			for (const auto& x : S) {
				if (x < lo) lo = x;
				if (x > hi) hi = x;
			}

			vec<compressed_value_type> buckets(hi - lo + 1, 0);
			for (const auto& x : S) {
				buckets[x - lo] = 1;
			}
			for (int v = 0; v < int(buckets.size()); v++) {
				if (buckets[v]) buckets[v] = compressed_value_type(sigma++);
			}

			for (int i = 0; i < int(std::ranges::size(S)); i++) {
				compressed_s[i] = buckets[s[i] - lo];
			}
		}

		return construct_raw(compressed_s, sigma);
	}

protected:
	template <integer_text String> void build(const String& S, index_t sigma) {
		N = int(std::ranges::size(S));
		build_sa(S, sigma);
		build_rank();
	}

private:
	template <integer_text String> void build_sa(const String& S, index_t sigma) {
		sa = vec<index_t>(N+1);
		assert(sigma >= 0);
		for (auto s : S) assert(0 <= index_t(s) && index_t(s) < sigma);
		vec<index_t> tmp(sigma + std::max(N, sigma));
		SuffixArrayBase::sais(N, std::ranges::begin(S), sa.data(), sigma, tmp.data());
	}

	template <std::random_access_iterator It> requires std::integral<std::iter_value_t<It>>
	static void sais(int N, It S, index_t* sa, int sigma, index_t* tmp) {
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
				sais(num_pieces, next_S, sa, next_sigma, tmp);

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
		rank = vec<index_t>(N+1, uninit);
		for (int i = 0; i <= N; i++) rank[sa[i]] = i;
	}
};

class SuffixArray : public SuffixArrayBase<SuffixArray> {};

template <typename Self> class SuffixArrayLCPBase : public SuffixArrayBase<Self> {
public:
	using index_t = SuffixArrayBase<Self>::index_t;
	// lcp[i] = lcp(sa[i], sa[i+1])
	vec<index_t> lcp;

protected:
	friend SuffixArrayBase<Self>;
	template <integer_text String> void build(const String& S, index_t sigma) {
		SuffixArrayBase<Self>::build(S, sigma);
		build_lcp(S);
	}

private:
	template <integer_text String> void build_lcp(const String& S) {
		auto s = std::ranges::begin(S);
		int N = this->N;
		const auto& sa = this->sa;
		const auto& rank = this->rank;
		assert(int(std::ranges::size(S)) == N);
		lcp = vec<index_t>(N);
		for (int i = 0, k = 0; i < N - 1; i++) {
			int j = sa[rank[i]-1];
			while (k < N - std::max(i, j) && s[i+k] == s[j+k]) k++;
			lcp[rank[i]-1] = k;
			if (k) --k;
		}
	}
};

class SuffixArrayLCP : public SuffixArrayLCPBase<SuffixArrayLCP> {};

template <typename Self> class SuffixArrayRMQBase : public SuffixArrayLCPBase<Self> {
public:
	using index_t = SuffixArrayLCPBase<Self>::index_t;
	RangeMinQuery<std::pair<index_t, index_t>> rmq;

	index_t get_lcp(index_t a, index_t b) const {
		if (a == b) return this->N-a;
		a = this->rank[a], b = this->rank[b];
		if (a > b) std::swap(a, b);
		return rmq.query(a, b-1).first;
	}

	// Get the split in the suffix tree, using half-open intervals
	// Returns len, idx
	std::pair<index_t, index_t> get_split(index_t l, index_t r) const {
		assert(r - l > 1);
		return rmq.query(l, r-2);
	}

protected:
	friend SuffixArrayBase<Self>;
	template <integer_text String> void build(const String& S, index_t sigma) {
		SuffixArrayLCPBase<Self>::build(S, sigma);
		build_rmq();
	}

private:
	void build_rmq() {
		int N = this->N;
		const auto& lcp = this->lcp;
		vec<std::pair<index_t, index_t>> lcp_idx(N);
		for (int i = 0; i < N; i++) {
			lcp_idx[i] = {lcp[i], i+1};
		}
		rmq = RangeMinQuery<std::pair<index_t, index_t>>(std::move(lcp_idx));
	}
};

class SuffixArrayRMQ : public SuffixArrayRMQBase<SuffixArrayRMQ> {};

class PrefixArrayRMQ : private SuffixArrayRMQ {
	PrefixArrayRMQ(const SuffixArrayRMQ& sa_) : SuffixArrayRMQ(sa_) {}
	PrefixArrayRMQ(SuffixArrayRMQ&& sa_) : SuffixArrayRMQ(std::move(sa_)) {}

public:
	PrefixArrayRMQ() {}
	template <integer_text String> static PrefixArrayRMQ construct_raw(const String& S, int sigma) {
		return PrefixArrayRMQ(SuffixArrayRMQ::construct_raw(vec<std::ranges::range_value_t<String>>(std::from_range, std::views::reverse(S)), sigma));
	}

	// TODO: Fill in other constructors

	int get_lcs(int a, int b) const {
		return SuffixArrayRMQ::get_lcp(N - a, N - b);
	}
};

} // namespace wala
