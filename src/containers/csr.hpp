#pragma once

#include <algorithm>
#include <ranges>
#include <span>
#include <utility>

#include "containers/vec.hpp"

namespace wala {

struct csr_index {
	vec<int> bounds;
	csr_index() = default;
	explicit csr_index(int num_rows) : bounds(num_rows + 1) {}
	explicit csr_index(vec<int> bounds_) : bounds(std::move(bounds_)) {}
	[[nodiscard]] std::ranges::iota_view<int, int> indices(int i) const { return std::views::iota(bounds[i], bounds[i+1]); }
	template <std::ranges::contiguous_range R> [[nodiscard]] auto slice(int i, R&& base) const {
		return std::span(base).subspan(bounds[i], bounds[i+1] - bounds[i]);
	}
	int num_rows() const { return bounds.empty() ? 0 : int(bounds.size()) - 1; }
	int num_entries() const { return bounds.empty() ? 0 : bounds.back(); }
};

template <typename T> struct csr : csr_index {
	vec<T> dat;
	csr() = default;
	template <typename... A> requires detail::uninitialized_constructor<T, A...>
	explicit csr(int num_rows, int num_entries, const A&... a) : csr_index(num_rows), dat(num_entries, a...) {}
	explicit csr(int num_rows, int num_entries, const T& v) requires std::copy_constructible<T> : csr_index(num_rows), dat(num_entries, v) {}
	explicit csr(csr_index idx, vec<T> dat_) : csr_index(std::move(idx)), dat(std::move(dat_)) {}
	std::span<T> operator [](int i) { return slice(i, dat); }
	std::span<const T> operator [](int i) const { return slice(i, dat); }
};

struct csr_index_builder {
	vec<int> bounds;
	csr_index_builder() = default;
	explicit csr_index_builder(int N) : bounds(N+1) {}
	void count(int k) { bounds[k+1]++; }
	csr_index finalize() && {
		for (int i = 1; i < int(bounds.size()); i++) {
			bounds[i] += bounds[i-1];
		}
		return csr_index(std::move(bounds));
	}
};

template <typename T> struct csr_builder {
	csr_index idx;
	vec<T> dat;
	csr_builder() = default;
	explicit csr_builder(csr_index idx_, vec<T>&& dat_buf = {}) : idx(std::move(idx_)), dat(std::move(dat_buf)) {
		int l = idx.num_entries();
		if (!idx.bounds.empty()) {
			std::shift_right(idx.bounds.begin(), idx.bounds.end(), 1);
			idx.bounds[0] = 0;
		}
		if (dat.size() != l) dat = vec<T>(l, uninit);
	}
	explicit csr_builder(csr_index_builder&& idx_builder, vec<T>&& dat_buf = {}) : idx(std::move(idx_builder.bounds)), dat(std::move(dat_buf)) {
		int l = 0;
		for (int i = 1; i < int(idx.bounds.size()); i++) {
			idx.bounds[i] = std::exchange(l, l + idx.bounds[i]);
		}
		if (dat.size() != l) dat = vec<T>(l, uninit);
	}
	[[nodiscard]] T& push(int k) { return dat[idx.bounds[k+1]++]; }
	[[nodiscard]] csr<T> finalize() && { return csr<T>(std::move(idx), std::move(dat)); }
};

} // namespace wala
