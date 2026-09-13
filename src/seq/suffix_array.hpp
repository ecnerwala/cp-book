#pragma once

/*
 * This is mostly inspired by https://golang.org/src/index/suffixarray/sais.go.
 */

#include <algorithm>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include <type_traits>

#include "ds/rmq.hpp"

namespace wala {


// Layered suffix array: SuffixArrayBase computes just sa, each further
// layer statically opts into one more derived structure. Use the leaf classes
// SuffixArray, SuffixArrayRank, SuffixArrayLCP, or SuffixArrayRMQ; the named
// constructors on each return that type.
template <typename Self> class SuffixArrayBase {
public:
	using index_t = int;
	int N;
	// sa[0] = N is the sentinel suffix.
	std::vector<index_t> sa;

	SuffixArrayBase() : N(0) {}

	template <typename String> static Self construct_raw(const String& S, index_t sigma) {
		Self res;
		res.build(S, sigma);
		return res;
	}

	// Pass a function which returns a value in [0, sigma)
	template <typename String, typename F> static Self map_and_construct(const String& S, const F& f, int sigma) {
		std::vector<decltype((f(S[0])))> mapped(int(std::size(S)));
		for (int i = 0; i < int(std::size(S)); i++) {
			mapped[i] = f(S[i]);
			assert(0 <= int(mapped[i]) && int(mapped[i]) < sigma);
		}
		return construct_raw(mapped, sigma);
	}

	// Sorts the elements of S and then runs suffix array. This takes O(N log N) time with no dependence on sigma.
	template <typename String> static Self sort_and_construct(const String& S) {
		using std::begin;
		using std::end;
		using value_type = typename std::iterator_traits<decltype(begin(S))>::value_type;
		using compressed_value_type = typename std::conditional<
			sizeof(value_type) < sizeof(index_t),
			value_type,
			index_t
		>::type;

		std::vector<compressed_value_type> compressed_s(int(std::size(S)));
		int sigma = 0;

		{
			std::vector<value_type> vals(begin(S), end(S));
			std::sort(vals.begin(), vals.end());
			vals.resize(unique(vals.begin(), vals.end()) - vals.begin());
			for (int i = 0; i < int(std::size(S)); i++) {
				compressed_s[i] = compressed_value_type(index_t(std::lower_bound(vals.begin(), vals.end(), S[i]) - vals.begin()));
			}
			sigma = int(vals.size());
		}

		return construct_raw(compressed_s, sigma);
	}

	// Shifts the elements so that sigma = max(S) - min(S) + 1
	template <typename String> static Self shift_and_construct(const String& S) {
		using std::begin;
		using std::end;
		using value_type = typename std::iterator_traits<decltype(begin(S))>::value_type;

		std::vector<value_type> compressed_s(int(std::size(S)));
		int sigma = 0;

		if (int(std::size(S)) > 0) {
			value_type lo = *begin(S), hi = *begin(S);
			for (const auto& x : S) {
				if (x < lo) lo = x;
				if (x > hi) hi = x;
			}

			for (int i = 0; i < int(std::size(S)); i++) {
				compressed_s[i] = value_type(S[i] - lo);
			}
			sigma = int(hi - lo + 1);
		}

		return construct_raw(compressed_s, sigma);
	}

	// Renumber/filter to only the used elements with bucket sorting. Still takes O(max(S) - min(S) + 1) memory/time,
	// but should be less memory than `shift_and_construct` when sigma ~ N and max(S) - min(S) + 1 > N.
	template <typename String> static Self bucket_and_construct(const String& S) {
		using std::begin;
		using std::end;
		using value_type = typename std::iterator_traits<decltype(begin(S))>::value_type;
		using compressed_value_type = typename std::conditional<
			sizeof(value_type) < sizeof(index_t),
			value_type,
			index_t
		>::type;

		std::vector<compressed_value_type> compressed_s(int(std::size(S)));
		int sigma = 0;

		if (int(std::size(S)) > 0) {
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

			for (int i = 0; i < int(std::size(S)); i++) {
				compressed_s[i] = buckets[S[i] - lo];
			}
		}

		return construct_raw(compressed_s, sigma);
	}

protected:
	template <typename String> void build(const String& S, index_t sigma) {
		N = int(std::size(S));
		assert(sigma >= 0);
		for (auto s : S) assert(0 <= index_t(s) && index_t(s) < sigma);
		sa = std::vector<index_t>(N+1);
		// Scratch for sais: the work array (N+1), then lms_pos (N/2+1).
		// The recursion (on at most N/2 pieces) fits in the work array.
		std::vector<index_t> tmp(N+1 + N/2+1);
		SuffixArrayBase::sais<String>(N, S, sa.data(), sigma, tmp.data());
	}

private:
	// Suffix array by induced sorting (SA-IS): computes sa[0..N] for S plus a sentinel.
	//
	// We classify each position by (own type, predecessor's type): A = L/L, B = L/S, D = S/S, C = S/L
	// (the LMS positions).
	// Inducing from an entry only ever produces something in the L pass for A and C entries and in the S
	// pass for B and D entries, so instead of scanning sa and skipping the rest, each pass reads exactly the
	// entries it processes from a work array W laid out as
	//   W = [A_0 C_0 A_1 C_1 ... A_{sigma-1} C_{sigma-1} | B_0 D_0 B_1 D_1 ... B_{sigma-1} D_{sigma-1}]
	// where X_c holds the class X entries of bucket c in sa order.
	// The L pass reads the first half front to back, and every A entry has been induced (from a smaller
	// A or C entry) by the time it's read; the S pass reads the second half back to front, likewise.
	// The L pass fills A and B front to back; the S pass fills D back to front, and writes the C entries
	// back to front into W[0..num_pieces), where they come out as the sorted list of LMS positions.
	// The LMS round uses sa itself as W; the final round has to write sa, so it uses tmp as W.
	template <typename String> static void sais(int N, const String& S, index_t* sa, int sigma, index_t* tmp) {
		if (N == 0) {
			sa[0] = 0;
			return;
		} else if (N == 1) {
			sa[0] = 1;
			sa[1] = 0;
			return;
		}

		index_t* const W_final = tmp;
		// LMS positions in decreasing order. This lives above W_final so that it survives the recursion.
		index_t* const lms_pos = tmp + (N+1);

		// hist[cls*sigma + c] counts the positions of class cls in bucket c, excluding position 0.
		// bkt[c] is the sa bucket pointer for the final round.
		// ptr[cls*sigma + c] is the write pointer into W for X_c, and ptr[4*sigma] a dead slot for position
		// 0, which nothing is induced from.
		// (Pointers rather than indices into W: the indexed stores that indices compile to measured 2x
		// slower on inputs with long same-bucket chains, presumably by defeating memory renaming.)
		std::vector<index_t> book(5*sigma);
		std::vector<index_t*> ptr(4*sigma + 1);
		index_t* const hist = book.data();
		index_t* const bkt = hist + 4*sigma;

		// Phase 1: classify, counting each class/bucket and recording the LMS positions.
		int num_pieces = 0;
		{
			index_t c0 = S[N-1], c1;
			bool isS0 = false, isS1;
			for (int i = N-2; i >= 0; i--) {
				c1 = c0, isS1 = isS0;
				c0 = S[i];
				isS0 = (c0 < c1) | ((c0 == c1) & isS1);
				bool lms = isS1 & !isS0;
				hist[(2*isS1 + (isS0 ^ isS1))*sigma + c1]++;
				lms_pos[num_pieces] = i+1;
				num_pieces += lms;
			}
		}

		// Phase 2: sort the LMS substrings, if there's more than one.
		if (num_pieces > 1) {
			induce<false>(N, S, sigma, hist, bkt, ptr.data(), sa, lms_pos, num_pieces, nullptr);
			index_t* const pieces = sa;

			// Compute the lengths of the pieces in preparation for equality
			// comparison, and store them in tmp[v/2]. We set the length of the
			// final piece to 0; it compares unequal to everything because of
			// the sentinel.
			tmp[lms_pos[0]>>1] = 0;
			for (int k = 1; k < num_pieces; k++) {
				int v = lms_pos[k];
				tmp[v>>1] = lms_pos[k-1] - v;
			}

			// Compute the alphabet, storing the result into tmp[v/2].
			int next_sigma = 0;
			{
				int prv_len = -1, prv_v = 0;
				for (int i = 0; i < num_pieces; i++) {
					int v = pieces[i];
					int len = tmp[v>>1];

					bool eq = prv_len == len;
					for (int a = 0; eq && a < len; ++a) {
						eq = S[v+a] == S[prv_v+a];
					}
					if (!eq) {
						next_sigma++;
						prv_len = len;
						prv_v = v;
					}

					tmp[v>>1] = next_sigma - 1;
				}
			}

			if (next_sigma == num_pieces) {
				memmove(sa+1, pieces, sizeof(*sa) * num_pieces);
			} else {
				// Pack the input to the recursion: the names in text order, at the top of sa.
				// The recursion's output and scratch stay below it.
				index_t* next_S = sa + N + 1 - num_pieces;
				for (int k = 0; k < num_pieces; k++) {
					next_S[num_pieces-1-k] = tmp[lms_pos[k]>>1];
				}

				sais<const index_t*>(num_pieces, next_S, sa, next_sigma, tmp);

				// Map the suffix array of the names back up to piece start points
				for (int i = 1; i <= num_pieces; i++) {
					sa[i] = lms_pos[num_pieces-1-sa[i]];
				}
			}
		} else if (num_pieces == 1) {
			sa[1] = lms_pos[0];
		}

		// Phase 3: induce everything from the sorted pieces, now in sa[1..num_pieces].
		induce<true>(N, S, sigma, hist, bkt, ptr.data(), W_final, sa+1, num_pieces, sa);
	}

	// One round of induced sorting from the given seeds (the LMS positions, sorted if FINAL).
	// If FINAL, also writes every suffix into sa in its final position.
	// Every entry v read induces u = v-1; only u = 0 induces nothing, which we route to the dead slot.
	template <bool FINAL, typename String> static void induce(
		int N, const String& S, int sigma,
		const index_t* hist, index_t* bkt, index_t** ptr,
		index_t* W, const index_t* seeds, int num_seeds, index_t* sa
	) {
		const index_t* const histA = hist;
		const index_t* const histB = hist + sigma;
		const index_t* const histD = hist + 2*sigma;
		const index_t* const histC = hist + 3*sigma;
		index_t** const ptrA = ptr;
		index_t** const ptrB = ptr + sigma;
		index_t** const ptrD = ptr + 2*sigma;
		index_t** const ptrC = ptr + 3*sigma;
		const int dead = 4*sigma;

		index_t* p = W;
		for (int c = 0; c < sigma; c++) {
			ptrA[c] = p, p += histA[c];
			ptrC[c] = p, p += histC[c];
		}
		index_t* const W_mid = p;
		for (int c = 0; c < sigma; c++) {
			ptrB[c] = p, p += histB[c] + histD[c];
			ptrD[c] = p;
		}
		assert(p == W + N-1);

		for (int i = 0; i < num_seeds; i++) {
			int v = seeds[i];
			*ptrC[index_t(S[v])]++ = v;
		}

		// L pass
		{
			if constexpr (FINAL) {
				int cur = 1;
				for (int c = 0; c < sigma; c++) {
					bkt[c] = cur;
					cur += histA[c] + histB[c] + histD[c] + histC[c] + (c == index_t(S[0]));
				}
				sa[0] = N;
			}
			ptr[dead] = W + N-1;

			auto push = [&](int u) {
				index_t c1 = S[u];
				index_t c0 = u ? S[u-1] : c1;
				bool predS = c0 < c1;
				if constexpr (FINAL) sa[bkt[c1]++] = u;
				int k = u ? predS*sigma + c1 : dead;
				*ptr[k]++ = u;
			};
			push(N-1);
			for (index_t* rd = W; rd < W_mid; rd++) push(*rd - 1);
		}

		// S pass
		{
			if constexpr (FINAL) {
				int cur = 1;
				for (int c = 0; c < sigma; c++) {
					cur += histA[c] + histB[c] + histD[c] + histC[c] + (c == index_t(S[0]));
					bkt[c] = cur;
				}
			}
			{
				index_t* q = W;
				for (int c = 0; c < sigma; c++) {
					q += histC[c];
					ptrC[c] = q;
				}
			}
			ptr[dead] = W + N;

			auto push = [&](int u) {
				index_t c1 = S[u];
				index_t c0 = u ? S[u-1] : c1+1;
				bool predL = c0 > c1;
				if constexpr (FINAL) sa[--bkt[c1]] = u;
				int k = u ? (2 + predL)*sigma + c1 : dead;
				*--ptr[k] = u;
			};
			for (index_t* rd = W + N-1; rd > W_mid; ) push(*--rd - 1);
		}
	}

};

class SuffixArray : public SuffixArrayBase<SuffixArray> {};

template <typename Self> class SuffixArrayRankBase : public SuffixArrayBase<Self> {
public:
	using index_t = typename SuffixArrayBase<Self>::index_t;
	// rank[sa[i]] = i
	std::vector<index_t> rank;

protected:
	friend SuffixArrayBase<Self>;
	template <typename String> void build(const String& S, index_t sigma) {
		SuffixArrayBase<Self>::build(S, sigma);
		build_rank();
	}

private:
	void build_rank() {
		int N = this->N;
		const auto& sa = this->sa;
		rank = std::vector<index_t>(N+1);
		for (int i = 0; i <= N; i++) rank[sa[i]] = i;
	}
};

class SuffixArrayRank : public SuffixArrayRankBase<SuffixArrayRank> {};

template <typename Self> class SuffixArrayLCPBase : public SuffixArrayRankBase<Self> {
public:
	using index_t = typename SuffixArrayRankBase<Self>::index_t;
	// lcp[i] = lcp(sa[i], sa[i+1])
	std::vector<index_t> lcp;

protected:
	friend SuffixArrayBase<Self>;
	template <typename String> void build(const String& S, index_t sigma) {
		SuffixArrayRankBase<Self>::build(S, sigma);
		build_lcp(S);
	}

private:
	template <typename String> void build_lcp(const String& S) {
		int N = this->N;
		const auto& sa = this->sa;
		const auto& rank = this->rank;
		assert(int(std::size(S)) == N);
		lcp = std::vector<index_t>(N);
		for (int i = 0, k = 0; i < N - 1; i++) {
			int j = sa[rank[i]-1];
			while (k < N - std::max(i, j) && S[i+k] == S[j+k]) k++;
			lcp[rank[i]-1] = k;
			if (k) --k;
		}
	}
};

class SuffixArrayLCP : public SuffixArrayLCPBase<SuffixArrayLCP> {};

template <typename Self> class SuffixArrayRMQBase : public SuffixArrayLCPBase<Self> {
public:
	using index_t = typename SuffixArrayLCPBase<Self>::index_t;
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
	template <typename String> void build(const String& S, index_t sigma) {
		SuffixArrayLCPBase<Self>::build(S, sigma);
		build_rmq();
	}

private:
	void build_rmq() {
		int N = this->N;
		const auto& lcp = this->lcp;
		std::vector<std::pair<index_t, index_t>> lcp_idx(N);
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
	template <typename String> static PrefixArrayRMQ construct_raw(const String& S, int sigma) {
		return PrefixArrayRMQ(SuffixArrayRMQ::construct_raw(String(S.rbegin(), S.rend()), sigma));
	}

	// TODO: Fill in other constructors

	int get_lcs(int a, int b) const {
		return SuffixArrayRMQ::get_lcp(N - a, N - b);
	}
};

} // namespace wala
