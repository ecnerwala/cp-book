#pragma once

#include <functional>
#include <vector>
#include <cassert>

template <typename T, class Compare = std::less<T>> class RangeMinQuery : private Compare {
	static const int BUCKET_SIZE = 16;
	static const int BUCKET_SIZE_LOG = 4;
	static_assert((BUCKET_SIZE & (BUCKET_SIZE-1)) == 0);
	static const int CACHE_LINE_ALIGNMENT = 64;
	const T* data = nullptr;
	int n = 0;
	std::vector<T> tables;
	T* pref_data = nullptr;
	T* suff_data = nullptr;
	T* sparse_table = nullptr;

private:
	static constexpr int num_buckets(int n) {
		return n >> BUCKET_SIZE_LOG;
	}
	static constexpr int num_levels(int n) {
		return num_buckets(n) ? 32 - __builtin_clz(num_buckets(n)) : 0;
	}
	static constexpr int table_size(int n) {
		return n + n + num_buckets(n) * num_levels(n);
	}
private:
	const T& min(const T& a, const T& b) const {
		return Compare::operator()(a, b) ? a : b;
	}
	void setmin(T& a, const T& b) const {
		if (Compare::operator()(b, a)) a = b;
	}

public:
	RangeMinQuery() {}
	template <typename Vec> explicit RangeMinQuery(const Vec& v, const Compare& comp_ = Compare()) : RangeMinQuery(v.data(), int(v.size()), comp_) {}
	RangeMinQuery(const T* data_, int n_, const Compare& comp_ = Compare())
		: Compare(comp_)
		, data(data_)
		, n(n_)
		, tables(table_size(n))
		, pref_data(tables.data())
		, suff_data(tables.data()+n)
		, sparse_table(tables.data()+n+n) {
		for (int i = 0; i < n; i++) {
			pref_data[i] = data[i];
			if (i & (BUCKET_SIZE-1)) {
				setmin(pref_data[i], pref_data[i-1]);
			}
		}
		for (int i = n-1; i >= 0; i--) {
			suff_data[i] = data[i];
			if (i+1 < n && ((i+1) & (BUCKET_SIZE-1))) {
				setmin(suff_data[i], suff_data[i+1]);
			}
		}
		for (int i = 0; i < num_buckets(n); i++) {
			sparse_table[i] = data[i * BUCKET_SIZE];
			for (int v = 1; v < BUCKET_SIZE; v++) {
				setmin(sparse_table[i], data[i * BUCKET_SIZE + v]);
			}
		}
		for (int l = 0; l+1 < num_levels(n); l++) {
			for (int i = 0; i + (1 << (l+1)) <= num_buckets(n); i++) {
				sparse_table[(l+1) * num_buckets(n) + i] = min(sparse_table[l * num_buckets(n) + i], sparse_table[l * num_buckets(n) + i + (1 << l)]);
			}
		}
	}

	T query(int l, int r) const {
		assert(l <= r);
		int bucket_l = (l >> BUCKET_SIZE_LOG);
		int bucket_r = (r >> BUCKET_SIZE_LOG);
		if (bucket_l == bucket_r) {
			T ans = data[l];
			for (int i = l+1; i <= r; i++) setmin(ans, data[i]);
			return ans;
		} else {
			T ans = min(suff_data[l], pref_data[r]);
			bucket_l++;
			if (bucket_l < bucket_r) {
				int level = (32 - __builtin_clz(bucket_r - bucket_l)) - 1;
				setmin(ans, sparse_table[level * num_buckets(n) + bucket_l]);
				setmin(ans, sparse_table[level * num_buckets(n) + bucket_r - (1 << level)]);
			}
			return ans;
		}
	}
};

template <typename T> using RangeMaxQuery = RangeMinQuery<T, std::greater<T>>;
