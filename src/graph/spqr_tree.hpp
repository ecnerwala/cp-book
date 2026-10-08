#pragma once

#include <algorithm>
#include <vector>
#include <array>
#include <span>
#include <utility>
#include <cassert>
#include <ranges>
#include <ostream>
#include <expected>
#include <type_traits>
#include <variant>
#include <memory>

namespace wala {

inline constexpr struct with_capacity_t {} with_capacity;
inline constexpr struct uninit_t {} uninit;

namespace detail {
template <typename T> T* allocate(int n) { return n ? std::allocator<T>().allocate(n) : nullptr; }
template <typename T> void deallocate(T* p, int n) { if (p) std::allocator<T>().deallocate(p, n); }
inline void check_index([[maybe_unused]] int i, [[maybe_unused]] int n) {
#ifdef _GLIBCXX_DEBUG
	assert(0 <= i && i < n);
#endif
}
} // namespace detail

template <typename T> struct bounded_vector;
template <typename T> struct bounded_stack;

// Heap array whose length is fixed at construction.
template <typename T> struct fixed_vector {
	using value_type = T;
	using size_type = int;
	using difference_type = std::ptrdiff_t;
	using iterator = T*;
	using const_iterator = const T*;

	T* base = nullptr;
	int sz = 0;

	fixed_vector() = default;
	~fixed_vector() { std::destroy_n(base, sz); detail::deallocate(base, sz); }

	friend void swap(fixed_vector& a, fixed_vector& b) noexcept { std::swap(a.base, b.base); std::swap(a.sz, b.sz); }
	fixed_vector(fixed_vector&& o) noexcept : fixed_vector() { swap(*this, o); }
	fixed_vector& operator= (fixed_vector&& o) noexcept { swap(*this, o); return *this; }
	fixed_vector(const fixed_vector&) = delete;
	fixed_vector& operator= (const fixed_vector&) = delete;

	explicit fixed_vector(int n) : base(detail::allocate<T>(n)), sz(n) { std::uninitialized_value_construct_n(base, n); }
	explicit fixed_vector(int n, const T& v) : base(detail::allocate<T>(n)), sz(n) { std::uninitialized_fill_n(base, n, v); }
	explicit fixed_vector(uninit_t, int n) : base(detail::allocate<T>(n)), sz(n) {
		static_assert(std::is_trivially_default_constructible_v<T>);
	}
	template <std::ranges::sized_range R>
	fixed_vector(std::from_range_t, R&& r) : base(detail::allocate<T>(int(std::ranges::size(r)))), sz(int(std::ranges::size(r))) {
		std::ranges::uninitialized_copy(r, std::span(base, sz));
	}

	[[nodiscard]] fixed_vector clone() const { return fixed_vector(std::from_range, *this); }

	int size() const { return sz; }
	[[nodiscard]] bool empty() const { return sz == 0; }
	T* data() { return base; }
	const T* data() const { return base; }
	T* begin() { return base; }
	T* end() { return base + sz; }
	const T* begin() const { return base; }
	const T* end() const { return base + sz; }
	T& operator[] (int i) { detail::check_index(i, sz); return base[i]; }
	const T& operator[] (int i) const { detail::check_index(i, sz); return base[i]; }
	T& front() { assert(sz > 0); return base[0]; }
	const T& front() const { assert(sz > 0); return base[0]; }
	T& back() { assert(sz > 0); return base[sz-1]; }
	const T& back() const { assert(sz > 0); return base[sz-1]; }

	// The result is full (size == capacity)
	[[nodiscard]] bounded_vector<T> into_bounded() &&;
	[[nodiscard]] bounded_stack<T> into_stack() &&;

	friend bool operator == (const fixed_vector& a, const fixed_vector& b) { return std::ranges::equal(a, b); }
};

// Vector whose capacity is fixed at construction; never reallocates.
template <typename T> struct bounded_vector {
	using value_type = T;
	using size_type = int;
	using difference_type = std::ptrdiff_t;
	using iterator = T*;
	using const_iterator = const T*;

	T* base = nullptr;
	int sz = 0;
	int cap = 0;

	bounded_vector() = default;
	~bounded_vector() { std::destroy_n(base, sz); detail::deallocate(base, cap); }

	friend void swap(bounded_vector& a, bounded_vector& b) noexcept { std::swap(a.base, b.base); std::swap(a.sz, b.sz); std::swap(a.cap, b.cap); }
	bounded_vector(bounded_vector&& o) noexcept : bounded_vector() { swap(*this, o); }
	bounded_vector& operator= (bounded_vector&& o) noexcept { swap(*this, o); return *this; }
	bounded_vector(const bounded_vector&) = delete;
	bounded_vector& operator= (const bounded_vector&) = delete;

	explicit bounded_vector(with_capacity_t, int n) : base(detail::allocate<T>(n)), sz(0), cap(n) {}
	template <std::ranges::sized_range R>
	bounded_vector(std::from_range_t, R&& r) : bounded_vector(with_capacity, int(std::ranges::size(r))) {
		std::ranges::uninitialized_copy(r, std::span(base, cap));
		sz = cap;
	}

	[[nodiscard]] bounded_vector clone() const {
		bounded_vector r(with_capacity, cap);
		std::uninitialized_copy_n(base, sz, r.base);
		r.sz = sz;
		return r;
	}

	int size() const { return sz; }
	int capacity() const { return cap; }
	[[nodiscard]] bool empty() const { return sz == 0; }
	[[nodiscard]] bool full() const { return sz == cap; }
	T* data() { return base; }
	const T* data() const { return base; }
	T* begin() { return base; }
	T* end() { return base + sz; }
	const T* begin() const { return base; }
	const T* end() const { return base + sz; }
	T& operator[] (int i) { detail::check_index(i, sz); return base[i]; }
	const T& operator[] (int i) const { detail::check_index(i, sz); return base[i]; }
	T& front() { assert(sz > 0); return base[0]; }
	const T& front() const { assert(sz > 0); return base[0]; }
	T& back() { assert(sz > 0); return base[sz-1]; }
	const T& back() const { assert(sz > 0); return base[sz-1]; }

	template <typename... Args>
	T& emplace_back(Args&&... args) { assert(sz < cap); return *std::construct_at(base + sz++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { assert(sz > 0); std::destroy_at(base + --sz); }
	void truncate(int n) { assert(0 <= n && n <= sz); std::destroy(base + n, base + sz); sz = n; }
	void truncate(T* nend) { truncate(int(nend - base)); }
	void clear() { truncate(0); }
	void grow_to(int n) { assert(sz <= n && n <= cap); std::uninitialized_value_construct(base + sz, base + n); sz = n; }
	void grow_to(int n, const T& v) { assert(sz <= n && n <= cap); std::uninitialized_fill(base + sz, base + n, v); sz = n; }
	void grow_to(uninit_t, int n) {
		static_assert(std::is_trivially_default_constructible_v<T>);
		assert(sz <= n && n <= cap); sz = n;
	}
	void assign(int n, const T& v) { clear(); grow_to(n, v); }

	[[nodiscard]] fixed_vector<T> into_full_fixed() && {
		assert(full());
		fixed_vector<T> r;
		r.base = std::exchange(base, nullptr);
		r.sz = std::exchange(sz, 0);
		cap = 0;
		return r;
	}
	[[nodiscard]] bounded_stack<T> into_stack() &&;

	friend bool operator == (const bounded_vector& a, const bounded_vector& b) { return std::ranges::equal(a, b); }
};

// bounded_vector with a pointer to the top instead of a size, for hot stacks.
template <typename T> struct bounded_stack {
	using value_type = T;
	using size_type = int;
	using difference_type = std::ptrdiff_t;
	using iterator = T*;
	using const_iterator = const T*;

	T* base = nullptr;
	T* top = nullptr;
	T* cap = nullptr;

	bounded_stack() = default;
	~bounded_stack() { std::destroy(base, top); detail::deallocate(base, int(cap - base)); }

	friend void swap(bounded_stack& a, bounded_stack& b) noexcept { std::swap(a.base, b.base); std::swap(a.top, b.top); std::swap(a.cap, b.cap); }
	bounded_stack(bounded_stack&& o) noexcept : bounded_stack() { swap(*this, o); }
	bounded_stack& operator= (bounded_stack&& o) noexcept { swap(*this, o); return *this; }
	bounded_stack(const bounded_stack&) = delete;
	bounded_stack& operator= (const bounded_stack&) = delete;

	explicit bounded_stack(with_capacity_t, int n) : base(detail::allocate<T>(n)), top(base), cap(base + n) {}
	template <std::ranges::sized_range R>
	bounded_stack(std::from_range_t, R&& r) : bounded_stack(with_capacity, int(std::ranges::size(r))) {
		std::ranges::uninitialized_copy(r, std::span(base, cap));
		top = cap;
	}

	[[nodiscard]] bounded_stack clone() const {
		bounded_stack r(with_capacity, capacity());
		r.top = std::uninitialized_copy(base, top, r.base);
		return r;
	}

	int size() const { return int(top - base); }
	int capacity() const { return int(cap - base); }
	[[nodiscard]] bool empty() const { return top == base; }
	[[nodiscard]] bool full() const { return top == cap; }
	T* data() { return base; }
	const T* data() const { return base; }
	T* begin() { return base; }
	T* end() { return top; }
	const T* begin() const { return base; }
	const T* end() const { return top; }
	T& operator[] (int i) { detail::check_index(i, size()); return base[i]; }
	const T& operator[] (int i) const { detail::check_index(i, size()); return base[i]; }
	T& front() { assert(top > base); return base[0]; }
	const T& front() const { assert(top > base); return base[0]; }
	T& back() { assert(top > base); return top[-1]; }
	const T& back() const { assert(top > base); return top[-1]; }

	template <typename... Args>
	T& emplace_back(Args&&... args) { assert(top < cap); return *std::construct_at(top++, std::forward<Args>(args)...); }
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { assert(top > base); std::destroy_at(--top); }
	void truncate(T* ntop) { assert(base <= ntop && ntop <= top); std::destroy(ntop, top); top = ntop; }
	void truncate(int n) { truncate(base + n); }
	void clear() { truncate(base); }
	void grow_to(int n) { assert(size() <= n && n <= capacity()); std::uninitialized_value_construct(top, base + n); top = base + n; }
	void grow_to(int n, const T& v) { assert(size() <= n && n <= capacity()); std::uninitialized_fill(top, base + n, v); top = base + n; }
	void grow_to(uninit_t, int n) {
		static_assert(std::is_trivially_default_constructible_v<T>);
		assert(size() <= n && n <= capacity()); top = base + n;
	}
	void assign(int n, const T& v) { clear(); grow_to(n, v); }

	[[nodiscard]] fixed_vector<T> into_full_fixed() && {
		assert(full());
		fixed_vector<T> r;
		r.sz = size();
		r.base = std::exchange(base, nullptr);
		top = cap = nullptr;
		return r;
	}
	[[nodiscard]] bounded_vector<T> into_vector() && {
		bounded_vector<T> r;
		r.sz = size();
		r.cap = capacity();
		r.base = std::exchange(base, nullptr);
		top = cap = nullptr;
		return r;
	}

	friend bool operator == (const bounded_stack& a, const bounded_stack& b) { return std::ranges::equal(a, b); }
};

template <typename T> bounded_vector<T> fixed_vector<T>::into_bounded() && {
	bounded_vector<T> r;
	r.sz = r.cap = std::exchange(sz, 0);
	r.base = std::exchange(base, nullptr);
	return r;
}
template <typename T> bounded_stack<T> fixed_vector<T>::into_stack() && {
	bounded_stack<T> r;
	r.base = std::exchange(base, nullptr);
	r.top = r.cap = r.base + std::exchange(sz, 0);
	return r;
}
template <typename T> bounded_stack<T> bounded_vector<T>::into_stack() && {
	bounded_stack<T> r;
	r.base = std::exchange(base, nullptr);
	r.top = r.base + std::exchange(sz, 0);
	r.cap = r.base + std::exchange(cap, 0);
	return r;
}

// std::vector equivalent: a bounded_vector that is replaced by a larger one when full.
template <typename T> struct growable_vector {
	using value_type = T;
	using size_type = int;
	using difference_type = std::ptrdiff_t;
	using iterator = T*;
	using const_iterator = const T*;

	bounded_vector<T> buf;

	growable_vector() = default;
	explicit growable_vector(with_capacity_t, int n) : buf(with_capacity, n) {}
	explicit growable_vector(int n) : buf(with_capacity, n) { buf.grow_to(n); }
	explicit growable_vector(int n, const T& v) : buf(with_capacity, n) { buf.grow_to(n, v); }
	explicit growable_vector(uninit_t, int n) : buf(with_capacity, n) { buf.grow_to(uninit, n); }
	template <std::ranges::sized_range R>
	growable_vector(std::from_range_t, R&& r) : buf(std::from_range, std::forward<R>(r)) {}
	explicit growable_vector(bounded_vector<T>&& b) : buf(std::move(b)) {}

	friend void swap(growable_vector& a, growable_vector& b) noexcept { swap(a.buf, b.buf); }
	[[nodiscard]] growable_vector clone() const { return growable_vector(buf.clone()); }

	int size() const { return buf.size(); }
	int capacity() const { return buf.capacity(); }
	[[nodiscard]] bool empty() const { return buf.empty(); }
	T* data() { return buf.data(); }
	const T* data() const { return buf.data(); }
	T* begin() { return buf.begin(); }
	T* end() { return buf.end(); }
	const T* begin() const { return buf.begin(); }
	const T* end() const { return buf.end(); }
	T& operator[] (int i) { return buf[i]; }
	const T& operator[] (int i) const { return buf[i]; }
	T& front() { return buf.front(); }
	const T& front() const { return buf.front(); }
	T& back() { return buf.back(); }
	const T& back() const { return buf.back(); }

	// Moves the elements into a fresh bounded_vector of capacity n.
	void reallocate(int n) {
		assert(n >= buf.size());
		bounded_vector<T> nbuf(with_capacity, n);
		nbuf.sz = buf.size();
		std::uninitialized_move_n(buf.base, buf.size(), nbuf.base);
		swap(buf, nbuf);
	}
	void reserve(int n) { if (n > buf.capacity()) reallocate(n); }
	void shrink_to_fit() { if (!buf.full()) reallocate(buf.size()); }
	int grown_capacity() const { return std::max(2 * buf.capacity(), 4); }

	template <typename... Args>
	T& emplace_back(Args&&... args) {
		if (!buf.full()) [[likely]] return buf.emplace_back(std::forward<Args>(args)...);
		// Construct before relocating: args may alias an element
		bounded_vector<T> nbuf(with_capacity, grown_capacity());
		T& r = *std::construct_at(nbuf.base + buf.size(), std::forward<Args>(args)...);
		std::uninitialized_move_n(buf.base, buf.size(), nbuf.base);
		nbuf.sz = buf.size() + 1;
		swap(buf, nbuf);
		return r;
	}
	void push_back(const T& v) { emplace_back(v); }
	void push_back(T&& v) { emplace_back(std::move(v)); }
	void pop_back() { buf.pop_back(); }
	void truncate(int n) { buf.truncate(n); }
	void truncate(T* nend) { buf.truncate(nend); }
	void clear() { buf.clear(); }
	void grow_to(int n) { reserve(n); buf.grow_to(n); }
	void grow_to(int n, const T& v) { reserve(n); buf.grow_to(n, v); }
	void grow_to(uninit_t, int n) { reserve(n); buf.grow_to(uninit, n); }
	void assign(int n, const T& v) { clear(); grow_to(n, v); }

	[[nodiscard]] bounded_vector<T> into_bounded() && { return std::move(buf); }
	[[nodiscard]] fixed_vector<T> into_fixed() && { shrink_to_fit(); return std::move(buf).into_full_fixed(); }

	friend bool operator == (const growable_vector& a, const growable_vector& b) { return a.buf == b.buf; }
};

struct csr_index {
	fixed_vector<int> bounds;
	std::ranges::iota_view<int, int> indices(int i) const { return std::views::iota(bounds[i], bounds[i+1]); }
	template <std::ranges::contiguous_range R> auto slice(int i, R&& base) const {
		return std::span(base).subspan(bounds[i], bounds[i+1] - bounds[i]);
	}
	int num_rows() const { return bounds.empty() ? 0 : int(bounds.size()) - 1; }
	int num_entries() const { return bounds.empty() ? 0 : bounds.back(); }
};

template <typename T> struct csr : csr_index {
	fixed_vector<T> dat;
	std::span<T> operator [](int i) { return slice(i, dat); }
	std::span<const T> operator [](int i) const { return slice(i, dat); }
};

struct csr_index_builder {
	fixed_vector<int> bounds;
	csr_index_builder() = default;
	explicit csr_index_builder(int N) : bounds(N+1) {}
	void count(int k) { bounds[k+1]++; }
	csr_index finalize() && {
		for (int i = 1; i < int(bounds.size()); i++) {
			bounds[i] += bounds[i-1];
		}
		return {std::move(bounds)};
	}
};

template <typename T> struct csr_builder {
	csr_index idx;
	fixed_vector<T> dat;
	csr_builder() = default;
	explicit csr_builder(csr_index idx_, fixed_vector<T>&& dat_buf = {}) : idx(std::move(idx_)), dat(std::move(dat_buf)) {
		int l = idx.num_entries();
		if (!idx.bounds.empty()) {
			std::shift_right(idx.bounds.begin(), idx.bounds.end(), 1);
			idx.bounds[0] = 0;
		}
		if (dat.size() != l) dat = fixed_vector<T>(uninit, l);
	}
	explicit csr_builder(csr_index_builder&& idx_builder, fixed_vector<T>&& dat_buf = {}) : idx{std::move(idx_builder.bounds)}, dat(std::move(dat_buf)) {
		int l = 0;
		for (int i = 1; i < int(idx.bounds.size()); i++) {
			idx.bounds[i] = std::exchange(l, l + idx.bounds[i]);
		}
		if (dat.size() != l) dat = fixed_vector<T>(uninit, l);
	}
	[[nodiscard]] T& push(int k) { return dat[idx.bounds[k+1]++]; }
	[[nodiscard]] csr<T> finalize() && { return { std::move(idx), std::move(dat) }; }
};

struct planar_spqr_tree;

struct spqr_tree {
	// The SPQR tree of a graph is a canonical/"maximal" decomposition of the graph by 2-vertex cuts.
	// The tree consists of nodes which are graphs of virtual edges (vedges), corresponding to nontrivial 2-vertex cuts.
	// Virtual edges are paired, and we can reassemble the graph by gluing nodes at their matching vedges (and removing the vedge).
	// Real edges are represented as special Q nodes which each contain exactly 1 real edge and exactly 1 vedge.
	//
	// Traditionally, the SPQR tree is defined for each biconnected component,
	// but we will embed the SPQR decompositions inside the block-cut tree to get a (rooted) decomposition of the entire graph.
	//
	// As such, we will have a tree of "items", which consist of SPQR nodes, real vertices, and a special "forest root" item:
	//  - Each vertex will be a child of the topmost node which contains it (or the forest root).
	//  - Each block will be a subtree of nodes rooted at a Q edge, which is the child of one of its vertices.
	//
	// Item types:
	//  F - forest - a root node corresponding to the whole forest.
	//  V - vertex - not really a node, just there because they're mixed into the tree a la block/cut tree
	//  Q - real edge - has exactly 1 vedge and 1 real edge
	//  I - bridge - has exactly 1 vedge connecting to a bridge Q node
	//  O - self-loop - has exactly 1 vedge connecting to a self-loop Q node
	//  S - series - a cycle of >= 3 vedges; note that any 2 vertices of the cycle form a cut
	//  P - parallel - a parallel group of >= 3 vedges with the same endpoints
	//  R - rigid - a 3-vertex-connected component
	//
	// Q nodes occur in 2 places: block roots and block leaves.
	// Block leaf Q's simply have no children.
	// Block root Q's have 2 children: their vedge, and their deeper vertex (unless it's a self-loop).
	//
	// Degenerate blocks:
	//  - a block consisting of a self-loop is a Q node connected to an O node.
	//  - a block consisting of a bridge is a Q node connected to an I node.
	//  - a block consisting of exactly 2 parallel edges is represented by 2 glued Q nodes.
	//
	// We have several id spaces:
	//  - items are in preorder
	//  - node_verts (nv's) are each node's vertices, given in node order then s-t order.
	//  - node_edges (ne's) are each node's vedges, given in node order then a s-t order.
	//  - node_adj is each node_vert's incident vedges, given as 2 lists per nv: left/rightwards edges each in reverse s-t order.
	//  - original verts and original edges can be converted to items as vert_item / edge_item
	//
	// Children of a node will be sorted in s-t order.
	// Specifically vertices are sorted, and edges are guaranteed to satisfy the strong "dominance" partial order:
	// if a.nvs[0] <= b.nvs[0] and a.nvs[1] <= b.nvs[1], then a <= b. (In practice, we'll sort by midpoint.)
	// Adjacency lists are sorted as "center-is-longest", which helps make laminar/bracket cases clean.
	//   (5->4) (5->3) (5->2) (5->1) *vertex 5* (5->9) (5->8) (5->7) (5->6)
	// More specifically, node_adj contains two lists per vertex: 2*nv+0 is leftwards and 2*nv+1 is rightwards.
	//
	// All id's are item indices unless clearly nv/ne id's.
	//
	// In general, there are 2 ways to use the SPQR tree: the rooted view and the unrooted view.
	//  - The rooted view uses par / ch walks, and either treats the tree as 1 top-down big decomposition, or walks in paths up/down the tree with LCA-like queries.
	//  - The unrooted view mostly uses nv/ne/nd lists and works locally within a node/sometimes jumps between them.

	enum class node_type : char {
		F = 'F', V = 'V', Q = 'Q', I = 'I', O = 'O', S = 'S', P = 'P', R = 'R'
	};
	friend std::ostream& operator<<(std::ostream& o, node_type t) { return o << char(t); }

	fixed_vector<int> vert_index;
	fixed_vector<int> edge_index;
	fixed_vector<bool> edge_flipped;

	fixed_vector<int> par;
	fixed_vector<int> subtree_end;
	fixed_vector<node_type> types;
	fixed_vector<int> orig_id;

	csr<int> ch;
	struct node_vert_t {
		int node;
		int vert;
	};
	fixed_vector<node_vert_t> node_verts;
	csr_index node_nvs;
	// The nv index of a vertex within its parent node
	fixed_vector<int> vert_par_nv;
	// TODO: Should we store a vert_nodes CSR?

	struct node_edge_t {
		int node;
		int twin_ne;
		// TODO: Should we store the twin node, the twin node type, and/or twin node type == Q?

		std::array<int, 2> nvs;
	};
	fixed_vector<node_edge_t> node_edges;
	csr_index node_nes;

	struct node_adj_t {
		int ne;
		int dest_nv;
	};
	csr<node_adj_t> node_adj;

	int size() const { return int(par.size()); }

	// vert_order and edge_order are (prefixes of) permutations of vertex / edge ids;
	// listed ids are visited first in the given order, then the rest in id order.
	// Roots are the first unvisited vertices, and DFS children are explored in edge order.
	// Use planar_spqr_tree::build to also compute the planar embeddings.
	static spqr_tree build(
		int NV,
		std::span<const std::array<int, 2>> edges,
		bool ternarize = false,
		std::span<const int> vert_order = {},
		std::span<const int> edge_order = {}
	) {
		return build_impl<false>(NV, edges, ternarize, vert_order, edge_order);
	}

protected:
	template <bool with_planarity>
	static std::conditional_t<with_planarity, planar_spqr_tree, spqr_tree> build_impl(
		int NV,
		std::span<const std::array<int, 2>> edges,
		bool ternarize,
		std::span<const int> vert_order,
		std::span<const int> edge_order
	);
};

struct planar_embedding {
	// We'll split our edges up into "quarter-edges", indexed according to
	//   4 * edge + 2 * side + dir, where side is v0 vs v1, and dir is cw vs ccw
	//
	//       1     2
	//    v0 ---e--- v1
	//       0     3
	//
	// We can think of a planar embedding as a collection of 3 involutions on quarter-edges:
	// * qe <-> qe ^ 1 maps quarter edges to their opposite side around the endpoint vertex.
	// * qe <-> qe ^ 3 maps quarter edges to their opposite side along the edge (around the face).
	// * qe <-> rot_adj[qe] maps quarter edges to their facing pair.
	// Walking around a vertex is alternating qe ^ 1 and rot_adj[qe], and walking around a face is qe ^ 3 and rot_adj[qe].
	//
	// Partial embeddings are represented with -1's in the rot_adj array.
	// NB: Helpers do not support -1's. It is up to the user to not access these entries!
	fixed_vector<int> rot_adj;
};

struct planar_spqr_tree : spqr_tree {
	fixed_vector<bool> node_planar;
	// Planarity adjacencies: ne_rot_adj is an involution of facing quarter-edges, indexed according to:
	// ne_rot_adj[4 * node_edge + 2 * side + dir]
	// Nonplanar nodes have all entries -1.
	planar_embedding ne_embedding;

	static planar_spqr_tree build(
		int NV,
		std::span<const std::array<int, 2>> edges,
		bool ternarize = false,
		std::span<const int> vert_order = {},
		std::span<const int> edge_order = {}
	) {
		return build_impl<true>(NV, edges, ternarize, vert_order, edge_order);
	}
};

// Phase 1: build a DFS skeleton with outedges sorted by lowval
struct lowval_storted_skeleton_t {
	bounded_vector<int> roots;
	struct key_t { int lowval; bool is_tree; bool is_type_1; };
	struct packed_key_t {
		int v;
		friend auto operator <=> (packed_key_t a, packed_key_t b) = default;
		[[nodiscard]] bool is_new_block() const { return v < 6; }
		[[nodiscard]] bool is_type_2() const { return v % 3 == 2; }
		[[nodiscard]] key_t unpack(int cur_depth) const {
			int lowval = v / 3 - 2; if (lowval < 0) lowval = cur_depth + ~lowval;
			int kind = v % 3;
			bool is_tree = kind != 1;
			bool is_type_1 = kind <= 1;
			return {lowval, is_tree, is_type_1 };
		}
	};
	struct outedge_t { int src, dest; int e_side; packed_key_t key; };
	csr<outedge_t> outedges;

	static lowval_storted_skeleton_t build(
		int NV,
		std::span<const std::array<int, 2>> edges,
		std::span<const int> vert_order,
		std::span<const int> edge_order
	) {
		// std::min is by reference, which breaks some optimizations
		auto min = [](auto a, auto b) { return a < b ? a : b; };

		int NE = int(edges.size());
		assert(int(vert_order.size()) <= NV);
		assert(int(edge_order.size()) <= NE);

		// Calls f(i) for i in order, then for the remaining i in [0, n) in increasing order.
		auto for_each_in_order = [] [[gnu::always_inline]] (int n, std::span<const int> order, auto f) -> void {
			for (int i : order) f(i);
			if (int(order.size()) == n) return;
			if (order.empty()) {
				for (int i = 0; i < n; i++) f(i);
			} else if (order.size() == 1) {
				for (int i = 0; i < n; i++) {
					if (i != order[0]) f(i);
				}
			} else {
				fixed_vector<bool> listed(n, false);
				for (int i : order) listed[i] = true;
				for (int i = 0; i < n; i++) {
					if (!listed[i]) f(i);
				}
			}
		};

		bounded_vector<int> roots(with_capacity, NV);
		csr<outedge_t> outedges;
		{
			fixed_vector<int> depth(NV, -1);
			// 1a: build a normal adjacency list for the initial lowval dfs
			struct edge_t { int dest; int e; };
			csr_index_builder adj_idx_builder(NV);
			for (auto [u, v] : edges) {
				adj_idx_builder.count(u);
				if (u != v) adj_idx_builder.count(v);
			}
			csr_builder<edge_t> adj_builder(std::move(adj_idx_builder));
			for_each_in_order(NE, edge_order, [&] [[gnu::always_inline]] (int e) -> void {
				auto [u, v] = edges[e];
				adj_builder.push(u) = {v, 2 * e + 0};
				if (u != v) adj_builder.push(v) = {u, 2 * e + 1};
			});
			auto adj = std::move(adj_builder).finalize();

			fixed_vector<outedge_t> all_outedges(NE);
			auto nxt_outedge = all_outedges.begin();
			// Return the 2 lowvals from this subtree
			struct dfs_stack_t {
				int cur;
				int prv_e;
				std::array<int, 2> lowvals;
				int ch_idx;
				int ch_end;
			};

			bounded_stack<dfs_stack_t> stk(with_capacity, NV);
			auto push_vert = [&] [[gnu::always_inline]] (int cur, int prv_e) -> void {
				int d = int(stk.size());
				depth[cur] = d;
				stk.push_back({cur, prv_e, {d, d}, adj.bounds[cur], adj.bounds[cur+1]});
			};
			auto finish_edge = [&] [[gnu::always_inline]] (bool is_tree, std::array<int, 2> n_lowvals) -> void {
				int d = int(stk.size()) - 1;
				auto& s = stk.back();
				int cur = s.cur;
				assert(s.ch_idx < s.ch_end);
				auto [nxt, e] = adj.dat[s.ch_idx];
				auto& lowvals = s.lowvals;
				s.ch_idx++;

				{
					// Extra bit is 0 for type-1 children, 1 for backedges, 2 for children with lowval2
					// Bridges have lowval -2 (kind 0), and components loops have lowval -1 (components are kind 0, loops are kind 1)
					// We group all backedges together/last to avoid breaking a straight-line graph embedding
					int lowval = n_lowvals[0];
					if (lowval >= d) lowval = ~(lowval - d);
					int kind = 2 * (n_lowvals[1] < d) + !is_tree;
					*nxt_outedge++ = {cur, nxt, e, packed_key_t{3 * (lowval + 2) + kind}};
				}

				// Keep the 2 distinct mins
				if (n_lowvals[0] < lowvals[0]) lowvals = {n_lowvals[0], min(n_lowvals[1], lowvals[0])};
				else lowvals[1] = min(lowvals[1], n_lowvals[0] == lowvals[0] ? n_lowvals[1] : n_lowvals[0]);
			};
			auto start_edge = [&] [[gnu::always_inline]] () -> void {
				int d = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);
				auto [nxt, e] = adj.dat[s.ch_idx];

				if ((e ^ 1) == s.prv_e || depth[nxt] > d) {
					// skip the edge
					s.ch_idx++; return;
				}

				bool is_tree = depth[nxt] == -1;
				if (is_tree) {
					push_vert(nxt, e);
				} else {
					finish_edge(false, {depth[nxt], d});
				}
			};
			auto pop_vert = [&] [[gnu::always_inline]] () -> std::array<int, 2> {
				auto lowvals = stk.back().lowvals;
				stk.pop_back();
				return lowvals;
			};
			for_each_in_order(NV, vert_order, [&] [[gnu::always_inline]] (int rt) -> void {
				if (depth[rt] == -1) {
					roots.push_back(rt);
					push_vert(rt, -1);
					while (true) {
						if (stk.back().ch_idx == stk.back().ch_end) {
							auto lowvals = pop_vert();
							if (stk.empty()) break;
							finish_edge(true, lowvals);
						} else {
							start_edge();
						}
					}
				}
			});

			assert(nxt_outedge == all_outedges.end());

			csr_index_builder by_key_idx_builder(3*NV+6);
			for (auto edge : all_outedges) by_key_idx_builder.count(edge.key.v);
			csr_builder<outedge_t> by_key_builder(std::move(by_key_idx_builder));
			for (auto edge : all_outedges) by_key_builder.push(edge.key.v) = edge;
			csr<outedge_t> by_key = std::move(by_key_builder).finalize();

			csr_index_builder by_src_idx_builder(NV);
			for (auto edge : by_key.dat) by_src_idx_builder.count(edge.src);
			csr_builder<outedge_t> by_src_builder(std::move(by_src_idx_builder), std::move(all_outedges));
			for (auto edge : by_key.dat) by_src_builder.push(edge.src) = edge;
			outedges = std::move(by_src_builder).finalize();
		}

		return {std::move(roots), std::move(outedges)};
	}
};

template <bool with_planarity>
std::conditional_t<with_planarity, planar_spqr_tree, spqr_tree> spqr_tree::build_impl(
	int NV,
	std::span<const std::array<int, 2>> edges,
	bool ternarize,
	std::span<const int> vert_order,
	std::span<const int> edge_order
) {
	// std::min is by reference, which breaks some optimizations
	auto setmin = [](auto& a, auto b) { if (b < a) a = b; };

	int NE = int(edges.size());
	assert(int(vert_order.size()) <= NV);
	assert(int(edge_order.size()) <= NE);

	auto [roots, outedges] = lowval_storted_skeleton_t::build(NV, edges, vert_order, edge_order);

	// Phase 2: do the big ear-decomposition-like walk

	// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
	// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

	// As we build, we will represent the children of our nodes/vertices as linked lists.
	constexpr int ROOT_ITEM = 0;
	auto vert_item = [&] [[gnu::always_inline]] (int v) -> int { return 1 + v; };
	auto edge_item = [&] [[gnu::always_inline]] (int e) -> int { return 1 + NV + e; };

	// Helpers for working with std::array<T, 2> - these compile to cmov's better than direct index access.

	// return arr[dir] == a, arr[!dir] == b
	auto set_sides = []<typename T>(bool dir, T a, T b) -> std::array<T, 2> {
		return dir ? std::array<T, 2>{b, a} : std::array<T, 2>{a, b};
	};
	auto get_side = []<typename T>(std::array<T, 2> a, bool dir) -> T {
		return dir ? a[1] : a[0];
	};

	struct item_list {
		// Items are actually 2 * item + planarity_flip (always 0 without planarity)
		std::array<int, 2> v{-1, -1};

		[[nodiscard]] bool empty() const { return v[0] < 0; }
	};
	bounded_vector<int> ch_nxt(with_capacity, 1 + NV + NE + NE);
	ch_nxt.grow_to(1 + NV + NE, -1);
	auto concat = [&] [[gnu::always_inline]] (item_list a, item_list b) -> item_list {
		if (b.empty()) return a;
		if (a.empty()) return b;
		ch_nxt[a.v[1] >> 1] = b.v[0] ^ (a.v[1] & 1);
		return {{a.v[0], b.v[1]}};
	};
	auto unit_list = [&] [[gnu::always_inline]] (int item) -> item_list {
		return {{item << 1, item << 1}};
	};

	bounded_vector<std::array<int, 2>> item_vs(with_capacity, 1 + NV + 2 * NE);
	item_vs.grow_to(1 + NV + NE, {-1, -1});
	bounded_vector<item_list> item_ch(with_capacity, 1 + NV + 2 * NE); item_ch.grow_to(1 + NV + NE, item_list{});
	bounded_vector<node_type> item_types(with_capacity, 1 + NV + 2 * NE);
	item_types.grow_to(1, node_type::F);
	item_types.grow_to(1 + NV, node_type::V);
	item_types.grow_to(1 + NV + NE, node_type::Q);

	// Quarter edges for planar embedding building.
	// Each vedge has 4 entries by 4 * vedge_id + 2 * source_vert + is_cw (is_cw is arbitrary)
	// vedges are identified with what item they cap, numbered by (item - 1 - NV)
	fixed_vector<int> quarter_edge_matches(with_planarity ? 8 * NE + 4 : 0, -1);
	struct nonplanarity_certficate_t {};
	bounded_vector<std::expected<std::array<int, 4>, nonplanarity_certficate_t>> node_planarity(with_capacity, with_planarity ? NE : 0);

	int tot_blocks = 0;
	int tot_self_loops = 0;

	{
		auto alloc_item = [&] [[gnu::always_inline]] (node_type type) -> int {
			int item = int(item_vs.size());
			item_vs.push_back({});
			item_ch.push_back({});
			item_types.push_back(type);
			ch_nxt.push_back(-1);
			if constexpr (with_planarity) node_planarity.emplace_back();
			return item;
		};

		// Declare these here: most of our code will be in terms of v_start / top_depth, so we'll want to read these out
		fixed_vector<int> stack_verts(NV);
		fixed_vector<bool> stack_dir(NV);

		auto make_vs = [&] [[gnu::always_inline]] (int v_start, int top_depth) -> std::array<int, 2> {
			return set_sides(stack_dir[top_depth], stack_verts[top_depth], v_start);
		};

		int nxt_edge_idx = 0; // Counts backedges only
		fixed_vector<int> first_occurrence(NV); // First backedge to this depth

		fixed_vector<int> edge_top_depths(with_planarity ? 2 * NE : 0, -1);

		struct tstack_planarity_side_t {
			// For each side, store pointers to the "linked lists" of the edges inside.
			// v[0] is the outer / longer edges and v[1] is the inner / shorter edges, matching the outside-in sort order.

			// bot_ends are the outer/innermost exposed pieces of the walk down the ear in the tree (they're connected to the bottommost/topmost vertices of the tree path)
			std::array<int, 2> bot_ends{-1, -1};
			struct top_t {
				int end = -1;
				int depth = -1;
			};
			// top_ends are the outer/innermost exposed backedges
			// depths should be increasing going inwards
			std::array<top_t, 2> tops{top_t{-1, -1}, top_t{-1, -1}};
		};
		struct tstack_planarity_t {
			// The convention is that sides[0].tops[0].depth == top_depth, i.e. at least one minimal return lives on side 0
			std::array<tstack_planarity_side_t, 2> sides;
		};
		struct tstack_nonplanarity_t {
			// TODO: What's the nonplanarity certificate look like?
		};
		using tstack_maybe_planarity_t = std::conditional_t<with_planarity, std::expected<tstack_planarity_t, tstack_nonplanarity_t>, std::monostate>;
		auto merge_planarity = [&] [[gnu::always_inline]] (tstack_maybe_planarity_t& a, const tstack_maybe_planarity_t& b) -> void {
			if constexpr (with_planarity) {
				if (!a) return;
				if (!b) { a = b; return; }
				for (int z = 0; z < 2; z++) {
					auto& as = a->sides[z];
					const auto& bs = b->sides[z];
					// If there's no bottom edges, then we must be an isolated vertex, so we can end early.
					if (bs.bot_ends[0] == -1) {
						// Do nothing
					} else if (as.bot_ends[0] == -1) {
						as = bs;
					} else {
						quarter_edge_matches[as.bot_ends[1]] = bs.bot_ends[0];
						quarter_edge_matches[bs.bot_ends[0]] = as.bot_ends[1];
						as.bot_ends[1] = bs.bot_ends[1];

						if (bs.tops[0].end == -1) {
							// Do nothing
						} else if (as.tops[0].end == -1) {
							as.tops = bs.tops;
						} else {
							// Caller must check that we're planar
							assert(as.tops[1].depth <= bs.tops[0].depth);
							quarter_edge_matches[as.tops[1].end] = bs.tops[0].end;
							quarter_edge_matches[bs.tops[0].end] = as.tops[1].end;
							as.tops[1] = bs.tops[1];
						}
					}
				}
			}
		};
		auto make_edge_planarity = [&] [[gnu::always_inline]] (int item, int top_depth, bool is_tree) -> tstack_maybe_planarity_t {
			if constexpr (with_planarity) {
				assert(item >= 1 + NV);
				int ve = item - (1 + NV);
				bool top_dir = stack_dir[top_depth];
				edge_top_depths[ve] = top_depth;
				tstack_planarity_t p;
				if (is_tree) {
					p.sides[0].bot_ends = {4 * ve + 2 * !top_dir + 0, 4 * ve + 2 * top_dir + 1};
					p.sides[1].bot_ends = {4 * ve + 2 * !top_dir + 1, 4 * ve + 2 * top_dir + 0};
				} else {
					p.sides[0].bot_ends = {4 * ve + 2 * !top_dir + 0, 4 * ve + 2 * !top_dir + 1};
					p.sides[0].tops = {{{4 * ve + 2 * top_dir + 1, top_depth}, {4 * ve + 2 * top_dir + 0, top_depth}}};
				}
				return p;
			} else {
				return {};
			}
		};
		struct tstack_t {
			int v_start = -1;
			int top_depth = -1;
			int first_idx = -1;
			std::array<item_list, 2> spans;
			[[no_unique_address]] tstack_maybe_planarity_t planarity;
		};
		bounded_stack<tstack_t> tstack(with_capacity, NV + NE);
		auto cur_tstack = [&] [[gnu::always_inline]] () -> tstack_t& { return tstack.end()[-1]; };
		auto nxt_tstack = [&] [[gnu::always_inline]] () -> tstack_t& { return tstack.end()[-2]; };

		auto push_tstack = [&] [[gnu::always_inline]] (int v_start, int top_depth, int item, tstack_maybe_planarity_t planarity) -> void {
			tstack.emplace_back(v_start, top_depth, nxt_edge_idx, set_sides(stack_dir[top_depth], unit_list(item), {}), planarity);
		};
		auto push_vert_tstack = [&] [[gnu::always_inline]] (int v, int top_depth) -> void {
			int item = vert_item(v);
			push_tstack(v, top_depth, item, {});
		};
		auto push_edge_tstack = [&] [[gnu::always_inline]] (int v_start, int top_depth, int e, bool is_tree) -> void {
			int item = edge_item(e);
			push_tstack(v_start, top_depth, item, make_edge_planarity(item, top_depth, is_tree));
		};
		auto flip_tstack_planarity = [&] [[gnu::always_inline]] (tstack_t& a) -> void {
			if constexpr (with_planarity) {
				a.spans[0].v[0] ^= 1;
				a.spans[0].v[1] ^= 1;
				a.spans[1].v[0] ^= 1;
				a.spans[1].v[1] ^= 1;
				if (a.planarity) {
					std::swap(a.planarity->sides[0], a.planarity->sides[1]);
				}
			}
		};
		auto merge_tstack_tops = [&] [[gnu::always_inline]] () -> void {
			tstack_t& a = nxt_tstack();
			const tstack_t& b = cur_tstack();
			setmin(a.top_depth, b.top_depth);
			a.spans[0] = concat(b.spans[0], a.spans[0]);
			a.spans[1] = concat(a.spans[1], b.spans[1]);
			if constexpr (with_planarity) {
				merge_planarity(a.planarity, b.planarity);
			}
			tstack.pop_back();
		};

		auto maybe_unwrap_nxt = [&] [[gnu::always_inline]] (node_type type, bool is_tree) -> int {
			tstack_t& t = nxt_tstack();

			if (type == node_type::R) return alloc_item(type);

			assert(type == node_type::P || type == node_type::S);

			// If we want to ternarize, never reuse.
			if (ternarize) return alloc_item(type);

			bool top_dir = stack_dir[t.top_depth];
			assert(get_side(t.spans, !top_dir).empty());
			int item = get_side(t.spans, top_dir).v[0] >> 1;
			assert(item == (get_side(t.spans, top_dir).v[1] >> 1));
			if (item_types[item] == type) {
				t.spans = set_sides(top_dir, item_ch[item], {});
				if constexpr (with_planarity) {
					// Unwrap the planarity data
					// We don't really need to maintain this at all because S/P nodes are known to be trivially planar
					// The current state is just make_edge_planarity(wrapped), which means that it has the right shape, just needs to be relabelled.
					assert(node_planarity[item - (1 + NV + NE)]);
					const auto& matches = *node_planarity[item - (1 + NV + NE)];
					assert(t.planarity);
					auto& p = *t.planarity;
					if (is_tree) {
						p.sides[0].bot_ends[0] = matches[2 * !top_dir + 1];
						p.sides[0].bot_ends[1] = matches[2 * top_dir + 0];
						p.sides[1].bot_ends[0] = matches[2 * !top_dir + 0];
						p.sides[1].bot_ends[1] = matches[2 * top_dir + 1];
					} else {
						p.sides[0].bot_ends[0] = matches[2 * !top_dir + 1];
						p.sides[0].bot_ends[1] = matches[2 * !top_dir + 0];
						p.sides[0].tops[0].end = matches[2 * top_dir + 0];
						p.sides[0].tops[1].end = matches[2 * top_dir + 1];
					}
				}
				return item;
			} else {
				return alloc_item(type);
			}
		};

		auto finish_tstack_top = [&] [[gnu::always_inline]] (int item, bool is_tree) -> void {
			tstack_t& t = cur_tstack();
			bool top_dir = stack_dir[t.top_depth];
			assert(get_side(t.spans, !top_dir).empty());

			if constexpr (with_planarity) {
				if (t.planarity) {
					const auto& p = *t.planarity;
					std::array<int, 4> matches{};
					if (is_tree) {
						matches[2 * !top_dir + 1] = p.sides[0].bot_ends[0];
						matches[2 * top_dir + 0] = p.sides[0].bot_ends[1];
						matches[2 * !top_dir + 0] = p.sides[1].bot_ends[0];
						matches[2 * top_dir + 1] = p.sides[1].bot_ends[1];
					} else {
						matches[2 * !top_dir + 1] = p.sides[0].bot_ends[0];
						matches[2 * !top_dir + 0] = p.sides[0].bot_ends[1];
						matches[2 * top_dir + 0] = p.sides[0].tops[0].end;
						matches[2 * top_dir + 1] = p.sides[0].tops[1].end;
					}
					node_planarity[item - (1 + NV + NE)] = matches;
				} else {
					assert(item_types[item] == node_type::R);
					node_planarity[item - (1 + NV + NE)] = std::unexpected(nonplanarity_certficate_t{});
				}
			}
			item_vs[item] = make_vs(t.v_start, t.top_depth);
			item_ch[item] = get_side(t.spans, top_dir);

			t.spans = set_sides(top_dir, unit_list(item), {});
			t.planarity = make_edge_planarity(item, t.top_depth, is_tree);
		};

		struct dfs_stack_t {
			bool has_vert_tstack;
			int ch_idx;
			int ch_end;
			int orig_tstack;
		};
		bounded_stack<dfs_stack_t> stk(with_capacity, NV);
		for (auto rt : roots) {
			auto push_vert = [&] [[gnu::always_inline]] (int cur) -> void {
				// stack_dir[cur_depth] must already be set for the lowval, so that we can push the vert tstack
				int cur_depth = int(stk.size());
				stack_verts[cur_depth] = cur;

				int lo = outedges.bounds[cur];
				int hi = outedges.bounds[cur+1];
				bool has_vert_tstack;
				{
					// Find the first same-BCC edge, and check it's type 2 (has lowval2), if so it's the ear tstack and we defer pushing ourselves.
					int first_edge = lo;
					while (first_edge < hi && outedges.dat[first_edge].key.is_new_block()) first_edge++;
					if (first_edge < hi && outedges.dat[first_edge].key.is_type_2()) {
						// Move first_edge to the beginning
						auto e = outedges.dat[first_edge];
						std::move_backward(outedges.dat.begin() + lo, outedges.dat.begin() + first_edge, outedges.dat.begin() + first_edge + 1);
						outedges.dat[lo] = e;
						has_vert_tstack = false;
					} else {
						push_vert_tstack(cur, cur_depth);
						has_vert_tstack = true;
					}
				}
				stk.push_back({has_vert_tstack, lo, hi, -1});
			};
			// return true means jump to start_edge, return false means jump to finish_edge
			auto start_edge = [&] [[gnu::always_inline]] () -> std::optional<int> {
				int cur_depth = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);
				auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
				auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

				// edge_dir convention: false is forwards, true is backwards.
				// That means that cur is on the edge_dir side and nxt is on the !edge_dir side.
				stack_dir[cur_depth] = (lowval >= cur_depth ? false : !stack_dir[lowval]);

				s.orig_tstack = int(tstack.size());
				if (is_tree) {
					stack_dir[cur_depth+1] = !stack_dir[std::min(lowval, cur_depth)];
					first_occurrence[cur_depth] = NE;
					return nxt;
				} else {
					return std::nullopt;
				}
			};
			auto finish_edge = [&] [[gnu::always_inline]] () -> void {
				int cur_depth = int(stk.size()) - 1;
				auto& s = stk.back();
				int cur = stack_verts[cur_depth];
				assert(s.ch_idx < s.ch_end);

				auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
				int e = e_side >> 1;
				s.ch_idx++;

				auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

				const int orig_tstack = s.orig_tstack;
				const bool edge_dir = stack_dir[cur_depth];

				if (lowval >= cur_depth) {
					// There's no planarity handling for this because it's just a Q node. I/O nodes also don't need any tracking.
					item_vs[edge_item(e)] = {cur, -1};
					tot_blocks++;
					if (is_tree) {
						// Bridges and components
						if (lowval == cur_depth + 1) {
							// tstack[tstack_size-1] is currently just smuggling out the child vertex, prepend the bridge component
							// This is just a shortcut for allocating a full I-type tstack
							int item = alloc_item(node_type::I);
							item_vs[item] = make_vs(nxt, cur_depth);
							item_ch[edge_item(e)] = concat(unit_list(item), tstack.back().spans[1]); tstack.pop_back();
						} else {
							// tstack[tstack_size-2] is the vertex and tstack[tstack_size-1] is the backedge
							auto backedge = tstack.back().spans[0]; tstack.pop_back();
							item_ch[edge_item(e)] = concat(backedge, tstack.back().spans[1]); tstack.pop_back();
						}
					} else {
						// self loops
						assert(nxt == cur);
						tot_self_loops++;
						int item = alloc_item(node_type::O);
						// Make sure the nxt is -1 as well
						item_vs[item] = {cur, -1};
						item_ch[edge_item(e)] = unit_list(item);
					}
					item_ch[vert_item(cur)] = concat(item_ch[vert_item(cur)], unit_list(edge_item(e)));
					return;
				}
				assert(lowval < cur_depth);

				item_vs[edge_item(e)] = make_vs(nxt, cur_depth);

				if (is_tree) {
					// The span lives on side edge_dir
					push_edge_tstack(nxt, cur_depth, e, true);
					while (nxt_tstack().top_depth >= cur_depth) {
						node_type type;
						if (nxt_tstack().top_depth > cur_depth) {
							// This is a vertex in the tstack, followed by either an S edge, possibly merged with other things
							if (tstack.end()[-3].top_depth < cur_depth) {
								// Not actually a good return, just stop
								break;
							}

							// Just backfill this for maybe_unwrap
							stack_dir[nxt_tstack().top_depth] = edge_dir;
							merge_tstack_tops();

							type = nxt_tstack().top_depth > cur_depth ? node_type::S : node_type::R;
						} else {
							assert(nxt_tstack().v_start == cur_tstack().v_start);
							// This will be a P node
							type = node_type::P;
						}
						int item = maybe_unwrap_nxt(type, type == node_type::S);
						merge_tstack_tops();
						if constexpr (with_planarity) {
							if (cur_tstack().planarity) {
								// Merge all backedges into the component
								for (auto& side : cur_tstack().planarity->sides) {
									assert(side.bot_ends[1] != -1);
									if (side.tops[1].end == -1) continue;
									assert(side.tops[0].depth == cur_depth);
									assert(side.tops[1].depth == cur_depth);
									quarter_edge_matches[side.bot_ends[1]] = side.tops[1].end;
									quarter_edge_matches[side.tops[1].end] = side.bot_ends[1];
									side.bot_ends[1] = side.tops[0].end;
									side.tops = {};
								}
							}
						}
						finish_tstack_top(item, true);
					}

					if (cur_tstack().first_idx > first_occurrence[cur_depth]) {
						if constexpr (with_planarity) {
							[&] [[gnu::always_inline]] () -> void {
								auto source = tstack.end() - 1;
								do {
									--source;
									if (!source->planarity) {
										// Set it somewhere so it can get copied down
										cur_tstack().planarity = source->planarity;
										return;
									}
								} while (source->first_idx > first_occurrence[cur_depth]);

								// From planarity's perspective, we can view each tstack as one of 2 shapes:
								// * tstack[i] can be a single "atom" branching off tstack[i].v_start. It can be:
								//   * A single backedge (type 1)
								//   * A subtree (possibly not biconnected) with at least 2 different-depth backedges on its outside (type 2)
								// * tstack[i] can be a "chunk". Chunks contain:
								//   * A core spanning from tstack[i].v_start (bot[0]) to tstack[i+1].v_start (bot[1]).
								//     * The core is a cyclic outer face: it has 2 *disjoint* paths from bot[0] to bot[1]
								//     * Each path can have backedges from its interior (not bot[0] or bot[1])
								//     * side[0]'s path has a backedge to tstack[i].lowval
								//   * Extra atoms on side 1, anchored at tstack[i].v_start
								//     * these must have lowval > tstack[i].lowval, or can have lowval == tstack[i].lowval and be type 2
								//     * these extra atoms are an entire suffix of v_start's: a chunk will always eat them all
								//   * Most of the time, we can treat the side[1] core and the atoms all as separate backedges from tstack[i].v_start.
								//     * The exception is when side[1].tops[1].depth == lowval: then it's forced to be a type 2 atom or part of the core, which matters.
								//       * TODO: Can we easily distinguish the 2 cases?
								// * tstack[i] can also be a tree vertex or a tree edge (trivial cases)
								//
								// Note that all atoms (including the chunk-extras) at one v_start must be sorted by (lowval, type).
								// However, a chunk can occur later (closer to the top) than its (lowval, type) sort at its v_start.

								// last_top == cur_tstack().top_depth
								int last_top = cur_depth;
								while (tstack.end() > source + 2) {
									if (nxt_tstack().top_depth > cur_depth) {
										// Vertex or tree edge, no conditions
									} else if (nxt_tstack().top_depth == cur_depth) {
										if (nxt_tstack().planarity->sides[1].tops[0].depth != -1) {
											// Double-sided to cur_depth, conflicts with cur_tstack()
											assert(last_top < cur_depth);
											// nxt_tstack() is a chunk and both backedges are on the core
											// K33 is:
											// * cur_tstack().tops[0]
											// * nxt_tstack().sides[0].tops[0].base
											// * nxt_tstack().sides[1].tops[0].base
											// + cur
											// + nxt_tstack().v_start
											// + cur_tstack().v_start
											//
											// nxt_tstack().v_start -> cur_tstack().tops[0] is the ear lowval loop
											// cur_tstack().v_start -> cur_tstack().tops[0] is just along cur
											// cur -> nxt_tstack().sides[*].tops[0].base is just the backedge
											// The rest is the outer face of nxt_tstack()
											cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
											return;
										}
										// We will put cur_depth on side 1 until the bottom
										flip_tstack_planarity(nxt_tstack());
									} else {
										if (nxt_tstack().planarity->sides[1].tops[0].depth != -1 && nxt_tstack().planarity->sides[1].tops[0].depth != cur_depth) {
											// Non-empty on both sides, conflicts with source
											// nxt_stack() is a chunk
											if (nxt_tstack().planarity->sides[1].tops[0].depth == nxt_tstack().top_depth) {
												// if it's core + type-2-atom
												// by the atom ordering, we're guaranteed source->cur isn't from v_start
												// * nxt_tstack().sides[0].tops[0].base
												// * nxt_tstack().sides[1].tops[0].base_fork
												// * source.base
												// + cur
												// + nxt_tstack().v_start
												// + nxt_tstack().top_depth
												//
												// cut nxt_tstack().core.side[1]
												// use nxt_tstack().sides[1].tops[0].prev to get from the fork to above cur down to cur
												//
												// otherwise it's double core
												// * cur
												// * nxt_tstack().sides[0].tops[0].base
												// * nxt_tstack().sides[1].tops[0].base
												// + cur_tstack().v_start
												// + nxt_tstack().v_start
												// + nxt_tstack().top_depth
												// (cut the lowval ear edge)
												//
												cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
											} else {
												// If it's an atom
												// by atom ordering, we're guaranteed source->cur isn't from v_start
												// * nxt_tstack().sides[0].tops[0].base
												// * nxt_tstack().sides[1].tops[0].end (go up/down to cur/top_depth)
												// * source.base
												// + cur
												// + nxt_tstack().v_start
												// + nxt_tstack().top_depth
												//
												// otherwise it's double core
												// * cur
												// * nxt_tstack().sides[0].tops[0].base
												// * nxt_tstack().sides[1].tops[0].base
												// + cur_tstack().v_start
												// + nxt_tstack().v_start
												// + nxt_tstack().sides[1].tops[0].end (side 0 gets there from above)
												// (cut the lowval ear edge)
												cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
											}
											return;
										}
										// Implicitly excludes -1
										if (nxt_tstack().planarity->sides[0].tops[1].depth > last_top) {
											assert(last_top < cur_depth);
											// 3 conflicting edges with nxt_tstack(), cur_tstack(), and source
											cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
											return;
										}
										last_top = nxt_tstack().top_depth;
									}
									merge_tstack_tops();
								}

								int t0 = nxt_tstack().planarity->sides[0].tops[1].depth;
								int t1 = nxt_tstack().planarity->sides[1].tops[1].depth;
								assert(t0 == cur_depth || t1 == cur_depth);
								if (std::min(t0, t1) > last_top) {
									assert(last_top < cur_depth);
									cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
									return;
								}
								if (t0 == cur_depth) {
									// We need to flip cur_tstack and nxt_tstack relative to each other.
									// Flip the one with worse top_depth.
									flip_tstack_planarity(cur_tstack().top_depth < nxt_tstack().top_depth ? nxt_tstack() : cur_tstack());
								}
								merge_tstack_tops();

								// Prune off finished cur-side things
								for (auto& side : cur_tstack().planarity->sides) {
									assert(side.bot_ends[1] != -1);
									while (side.tops[1].depth == cur_depth) {
										{
											// Link these to bot_ends[1]
											quarter_edge_matches[side.bot_ends[1]] = side.tops[1].end;
											quarter_edge_matches[side.tops[1].end] = side.bot_ends[1];
											side.bot_ends[1] = side.tops[1].end ^ 1;
										}
										side.tops[1].end = std::exchange(quarter_edge_matches[side.bot_ends[1]], -1);
										if (side.tops[1].end != -1) {
											quarter_edge_matches[side.tops[1].end] = -1;
											side.tops[1].depth = edge_top_depths[side.tops[1].end >> 2];
										} else {
											side.tops = {};
										}
									}
								}
							}();
						}
						while (cur_tstack().first_idx > first_occurrence[cur_depth]) {
							merge_tstack_tops();
						}
					}

					if (is_type_1) assert(s.has_vert_tstack);
					if (s.has_vert_tstack) {
						// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
						assert(int(tstack.size()) >= orig_tstack + 3);

						if (!is_type_1) {
							if constexpr (with_planarity) {
								[&] [[gnu::always_inline]] () -> void {
									// The lowval side should be side 1, everything else goes on side 0.
									// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
									// but is guaranteed to have *something* > lowval by non-type-1-ness
									auto& t = tstack[orig_tstack + 2];
									if (!t.planarity) {
										cur_tstack().planarity = t.planarity;
										return;
									}
									assert(t.planarity->sides[0].tops[0].depth == t.top_depth);
									assert(t.planarity->sides[0].tops[1].depth != -1);
									if (t.planarity->sides[0].tops[1].depth == lowval) {
										flip_tstack_planarity(t);
									} else if (t.planarity->sides[1].tops[1].depth != -1 && t.planarity->sides[1].tops[1].depth != lowval) {
										cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
										return;
									}
									assert(t.planarity->sides[0].tops[1].depth > lowval);
									int last_top = t.planarity->sides[0].tops[1].depth;
									for (int i = orig_tstack + 3; i < int(tstack.size()); i++) {
										if (!tstack[i].planarity) {
											cur_tstack().planarity = tstack[i].planarity;
											return;
										}
										if (tstack[i].top_depth == lowval) {
											flip_tstack_planarity(tstack[i]);
										}
										if (tstack[i].planarity->sides[1].tops[1].depth != -1 && tstack[i].planarity->sides[1].tops[1].depth != lowval) {
											cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
											return;
										}
										int next_top = tstack[i].planarity->sides[0].tops[0].depth;
										if (next_top != -1) {
											if (last_top > next_top) {
												cur_tstack().planarity = std::unexpected(tstack_nonplanarity_t{});
												return;
											}
											last_top = tstack[i].planarity->sides[0].tops[1].depth;
										}
									}
								}();
							}
							while (int(tstack.size()) > orig_tstack + 3) {
								merge_tstack_tops();
							}
						}

						assert(int(tstack.size()) == orig_tstack + 3);
						int item;
						if (is_type_1) {
							item = maybe_unwrap_nxt(cur_tstack().top_depth == cur_depth ? node_type::S : node_type::R, false);
						} else {
							// Just for the type checker
							item = -1;
						}
						// Merge with the backedge
						merge_tstack_tops();
						// Merge with the vertex
						merge_tstack_tops();

						cur_tstack().v_start = cur;
						assert(cur_tstack().top_depth == lowval);

						// Fold everything to the correct side now that we're leaving the child.
						// The entire subtree should go to the !edge_dir side.
						cur_tstack().spans = set_sides(!edge_dir, concat(cur_tstack().spans[0], cur_tstack().spans[1]), {});

						if constexpr (with_planarity) {
							if (cur_tstack().planarity) {
								// precondition: side 1 should be the lowval only side
								auto& sides = cur_tstack().planarity->sides;
								auto& s0 = sides[0];
								auto& s1 = sides[1];
								quarter_edge_matches[s0.bot_ends[0]] = s1.bot_ends[0];
								quarter_edge_matches[s1.bot_ends[0]] = s0.bot_ends[0];
								s0.bot_ends[0] = s1.bot_ends[1];
								if (s1.tops[0].end != -1) {
									// Caller must have checked that we're planar
									assert(s1.tops[1].depth == lowval);

									// This is always true
									assert(s1.tops[0].depth == lowval);
									quarter_edge_matches[s0.tops[0].end] = s1.tops[0].end;
									quarter_edge_matches[s1.tops[0].end] = s0.tops[0].end;
									s0.tops[0].end = s1.tops[1].end;
									// Already true since the backedge was on side 0
									assert(s0.tops[0].depth == lowval);
								}
								s1 = tstack_planarity_side_t{};
							}
						}

						if (is_type_1) {
							finish_tstack_top(item, false);
						}
					}
				} else {
					assert(is_type_1);
					// The span lives on side !edge_dir
					push_edge_tstack(cur, lowval, e, false);
					setmin(first_occurrence[lowval], nxt_edge_idx++);
				}

				assert(int(tstack.size()) >= orig_tstack + 1);

				// If is_type_1, the last entry on the tstack is either the vert_tstack, or the previous child as a unit
				if (is_type_1 && nxt_tstack().top_depth == lowval) {
					// This will be a P node
					int item = maybe_unwrap_nxt(node_type::P, false);
					merge_tstack_tops();
					finish_tstack_top(item, false);
				}

				if (!s.has_vert_tstack) {
					// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
					push_vert_tstack(cur, cur_depth);
					s.has_vert_tstack = true;
					assert(!is_type_1);
				}
			};
			auto pop_vert = [&] [[gnu::always_inline]] () -> void {
				auto& s = stk.back();
				assert(s.ch_idx == s.ch_end);
				assert(s.has_vert_tstack);
				stk.pop_back();
			};

			// Set something arbitrary, this is the normal convention for block-roots
			stack_dir[0] = true;
			push_vert(rt);
			while (true) {
				if (stk.back().ch_idx == stk.back().ch_end) {
					pop_vert();
					if (stk.empty()) break;
					finish_edge();
				} else if (std::optional<int> nxt = start_edge(); nxt) {
					push_vert(*nxt);
				} else {
					finish_edge();
				}
			}
			item_ch[ROOT_ITEM] = concat(item_ch[ROOT_ITEM], tstack.back().spans[1]); tstack.pop_back();
		}
	}

	// Phase 3: relabel the full tree in preorder
	int tot_items = int(item_types.size());
	{
		fixed_vector<int> vert_index(NV, -1);
		fixed_vector<int> edge_index(NE, -1);
		fixed_vector<bool> edge_flipped(NE, false);

		fixed_vector<int> par(tot_items, -1);
		fixed_vector<int> subtree_end(tot_items, -1);
		fixed_vector<node_type> types(tot_items, node_type::F);
		fixed_vector<int> orig_id(tot_items, -1);

		csr<int> ch;
		ch.bounds = fixed_vector<int>(tot_items + 1, 0);
		ch.dat = fixed_vector<int>(tot_items - 1);

		// Each node is a child, and additionally most non-block node has 2 cap verts; blocks have 1, and O nodes have 1
		int tot_node_verts = NV + (tot_items - 1 - NV) * 2 - tot_blocks - tot_self_loops;
		fixed_vector<node_vert_t> node_verts(tot_node_verts);
		csr_index node_nvs; node_nvs.bounds = fixed_vector<int>(tot_items + 1);
		fixed_vector<int> vert_par_nv(tot_items, -1);

		int tot_node_edges = (tot_items - 1 - NV - tot_blocks) * 2;
		fixed_vector<node_edge_t> node_edges(tot_node_edges);
		csr_index node_nes; node_nes.bounds = fixed_vector<int>(tot_items + 1);

		csr<node_adj_t> node_adj;
		node_adj.bounds = fixed_vector<int>(tot_node_verts * 2 + 1);
		node_adj.dat = fixed_vector<node_adj_t>(tot_node_edges * 2);

		fixed_vector<bool> node_planar(with_planarity ? tot_items : 0);
		fixed_vector<int> ne_rot_adj(with_planarity ? 4 * tot_node_edges : 0, -1);

		fixed_vector<int> vert_pos_buf(NV, -1);
		bounded_vector<int> cnts_buf(with_capacity, 2 * NV);
		struct ch_buf_t {
			int loc;
			int item_id;
		};
		bounded_vector<ch_buf_t> ch_buf(with_capacity, tot_items);
		fixed_vector<int> rot_edge_ne(with_planarity ? 2 * NE + 1 : 0);

		int nxt_unassigned_idx = 0;

		struct dfs_stack_t {
			int cur_idx;
			int ch_idx;
			int ch_end;
			int cur_nv;
			int cur_ne;
		};
		bounded_stack<dfs_stack_t> stk(with_capacity, tot_items);
		auto push_item = [&] [[gnu::always_inline]] (int cur_item) -> void {
			int cur_idx = nxt_unassigned_idx++;
			node_type cur_type = types[cur_idx] = item_types[cur_item];
			bool planar = true;
			if (cur_type == node_type::F) {
				assert(cur_item == 0);
			} else if (cur_type == node_type::V) {
				assert(1 <= cur_item && cur_item < 1 + NV);
				int orig_vert = cur_item - 1;
				orig_id[cur_idx] = orig_vert;
				vert_index[orig_vert] = cur_idx;
			} else if (cur_type == node_type::Q) {
				assert(1 + NV <= cur_item && cur_item < 1 + NV + NE);
				int orig_edge = cur_item - 1 - NV;
				orig_id[cur_idx] = orig_edge;
				edge_index[orig_edge] = cur_idx;
				assert(item_vs[cur_item][0] != -1);
				edge_flipped[orig_edge] = item_vs[cur_item][0] != edges[orig_edge][0];
			} else {
				assert(1 + NV + NE <= cur_item);
				if constexpr (with_planarity) {
					if (cur_type == node_type::O || cur_type == node_type::I) {
						// No planarity data was set up
					} else if (cur_type == node_type::S || cur_type == node_type::P || cur_type == node_type::R) {
						const auto& p = node_planarity[cur_item - (1 + NV + NE)];
						if (p) {
							// Make sure this runs before our planarity_flip checks
							for (int s = 0; s < 4; s++) {
								int a = 8 * NE + s, b = (*p)[s];
								quarter_edge_matches[a] = b;
								quarter_edge_matches[b] = a;
							}
						} else {
							// TODO: Any certificate stuff
							planar = false;
						}
					} else assert(false);
				}
			}
			if constexpr (with_planarity) node_planar[cur_idx] = planar;

			// HACK: Fill ch and vert_items in with orig items / orig verts for now,
			// because we don't have the final item id's yet.
			int ch_st = ch.bounds[cur_idx];
			int ch_en = ch_st;
			int nv_st = node_nvs.bounds[cur_idx];
			int nv_en = nv_st;
			int n_edges = 0;
			if (item_vs[cur_item][0] != -1) {
				node_verts[nv_en++] = {cur_idx, item_vs[cur_item][0]};
			}
			if (!item_ch[cur_item].empty()) {
				bool planarity_flip = item_ch[cur_item].v[0] & 1;
				for (int ch_item = item_ch[cur_item].v[0] >> 1; true; planarity_flip ^= (ch_nxt[ch_item] & 1), ch_item = ch_nxt[ch_item] >> 1) {
					ch.dat[ch_en++] = ch_item;
					assert(ch_item >= 1);
					if (ch_item < 1 + NV) {
						node_verts[nv_en++] = {cur_idx, ch_item - 1};
					} else {
						if constexpr (with_planarity) {
							if (cur_type != node_type::R) {
								assert(!planarity_flip);
							} else {
								// Fix the planarity direction right here: reverse quarter_edge_matches upfront;
								// this breaks the involution property, but from here on we'll never read the low bits anyways.
								int ve = ch_item - (1 + NV);
								if (planarity_flip) {
									std::swap(quarter_edge_matches[4 * ve + 0], quarter_edge_matches[4 * ve + 1]);
									std::swap(quarter_edge_matches[4 * ve + 2], quarter_edge_matches[4 * ve + 3]);
								}
							}
						}
						n_edges++;
					}
					if (ch_item == (item_ch[cur_item].v[1] >> 1)) {
						assert(ch_nxt[ch_item] == -1);
						break;
					}
				}
				planarity_flip ^= item_ch[cur_item].v[1] & 1;
				assert(!planarity_flip);
			}
			if (item_vs[cur_item][1] != -1) {
				node_verts[nv_en++] = {cur_idx, item_vs[cur_item][1]};
			}
			ch.bounds[cur_idx+1] = ch_en;
			node_nvs.bounds[cur_idx+1] = nv_en;

			int n_verts = nv_en - nv_st;

			bool is_node = cur_type != node_type::F && cur_type != node_type::V;
			bool has_cap = is_node && !(cur_type == node_type::Q && ch_en - ch_st > 0);

			if (!is_node) n_edges = 0;
			if (has_cap) n_edges++;

			int ne_st = node_nes.bounds[cur_idx];
			int ne_en = node_nes.bounds[cur_idx+1] = ne_st + n_edges;

			auto set_ne = [&] [[gnu::always_inline]] (int ne, std::array<int, 2> nvs, std::array<int, 2> nds, std::array<int, 4> rot_adjs) -> void {
				node_edges[ne].node = cur_idx;
				node_edges[ne].nvs = nvs;
				node_adj.dat[nds[0]] = {ne, nvs[1]};
				node_adj.dat[nds[1]] = {ne, nvs[0]};
				if constexpr (with_planarity) {
					for (int z = 0; z < 4; z++) ne_rot_adj[4 * ne + z] = rot_adjs[z];
				}
			};
			if (cur_type == node_type::F) {
				// Just set node_adj bounds and we're good
				for (int i = 2 * nv_st+1; i <= 2 * nv_en; i++) {
					node_adj.bounds[i] = 2 * ne_st;
				}
			} else if (cur_type == node_type::V) {
				// Nothing to do
			} else if (n_verts == 1) {
				assert(cur_type == node_type::Q || cur_type == node_type::O);
				assert(n_edges == 1);
				node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 1 * n_edges;
				node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 2 * n_edges;
				set_ne(ne_st, {nv_st, nv_st}, {2 * ne_st + 1, 2 * ne_st}, {4 * ne_st + 3, 4 * ne_st + 2, 4 * ne_st + 1, 4 * ne_st + 0});
			} else if (cur_type == node_type::Q || cur_type == node_type::I) {
				assert(n_verts == 2);
				assert(n_edges == 1);
				node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 0 * n_edges;
				node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 1 * n_edges;
				node_adj.bounds[2 * nv_st + 3] = 2 * ne_st + 2 * n_edges;
				node_adj.bounds[2 * nv_st + 4] = 2 * ne_st + 2 * n_edges;
				set_ne(ne_st, {nv_st, nv_st + 1}, {2 * ne_st, 2 * ne_st + 1}, {4 * ne_st + 1, 4 * ne_st + 0, 4 * ne_st + 3, 4 * ne_st + 2});
			} else if (cur_type == node_type::P) {
				// Special case: tiebreak the parallel edges so they're reversed
				assert(n_verts == 2);
				assert(n_edges >= 3);
				node_adj.bounds[2 * nv_st + 1] = 2 * ne_st + 0 * n_edges;
				node_adj.bounds[2 * nv_st + 2] = 2 * ne_st + 1 * n_edges;
				node_adj.bounds[2 * nv_st + 3] = 2 * ne_st + 2 * n_edges;
				node_adj.bounds[2 * nv_st + 4] = 2 * ne_st + 2 * n_edges;
				for (int ne = ne_st; ne < ne_en; ne++) {
					int ne_prv = (ne == ne_st ? ne_en : ne) - 1;
					int ne_nxt = (ne+1 == ne_en ? ne_st : ne+1);
					std::array<int, 4> rot_adjs{4 * ne_prv + 1, 4 * ne_nxt + 0, 4 * ne_nxt + 3, 4 * ne_prv + 2};
					set_ne(ne, {nv_st, nv_st + 1}, {2 * ne_st + (ne - ne_st), 2 * ne_en - 1 - (ne - ne_st)}, rot_adjs);
				}
			} else if (cur_type == node_type::S) {
				assert(n_verts == n_edges);
				assert(n_verts >= 3);
				for (int i = 2 * nv_st + 1; i <= 2 * nv_en; i++) {
					node_adj.bounds[i] = i + 2 * (ne_st - nv_st);
				}
				// Fix bounds for the cap
				node_adj.bounds[2 * nv_st + 1]--;
				node_adj.bounds[2 * nv_en - 1]++;
				set_ne(ne_st, {nv_st, nv_en - 1}, {2 * ne_st, 2 * ne_en - 1}, {4 * (ne_st+1) + 1, 4 * (ne_st+1) + 0, 4 * (ne_en-1) + 3, 4 * (ne_en-1) + 2});
				for (int i = 1; i < n_edges; i++) {
					int ne = ne_st + i;
					std::array<int, 4> rot_adjs{4 * (ne-1) + 3, 4 * (ne-1) + 2, 4 * (ne+1) + 1, 4 * (ne+1) + 0};
					if (ne-1 == ne_st) { rot_adjs[0] = 4 * ne_st + 1, rot_adjs[1] = 4 * ne_st + 0; }
					if (ne+1 == ne_en) { rot_adjs[2] = 4 * ne_st + 3, rot_adjs[3] = 4 * ne_st + 2; }
					set_ne(ne, {nv_st + i - 1, nv_st + i}, {2 * ne - 1, 2 * ne}, rot_adjs);
				}
			} else if (cur_type == node_type::R) {
				// Bucketsort the children by the midpoint
				for (int nv = nv_st; nv < nv_en; nv++) {
					vert_pos_buf[node_verts[nv].vert] = nv;
				}
				cnts_buf.assign(n_verts * 2 - 1, 0);
				ch_buf.clear();

				assert(has_cap);

				// Cap node_adj bounds
				node_adj.bounds[2 * nv_st + 2]++;
				node_adj.bounds[2 * nv_en - 1]++;

				for (int i = ch_st; i < ch_en; i++) {
					int item = ch.dat[i];
					assert(item >= 1);
					std::array<int, 2> nvs;
					if (item < 1 + NV) {
						nvs = {vert_pos_buf[item-1], vert_pos_buf[item-1]};
					} else {
						nvs = {vert_pos_buf[item_vs[item][0]], vert_pos_buf[item_vs[item][1]]};
						assert(nvs[0] < nvs[1]);
						node_adj.bounds[2 * nvs[0] + 2]++;
						node_adj.bounds[2 * nvs[1] + 1]++;
					}
					int loc = (nvs[0] - nv_st) + (nvs[1] - nv_st);
					ch_buf.emplace_back(loc, item);
					cnts_buf[loc]++;
				}
				int offset = ch_st;
				for (auto& cnt : cnts_buf) {
					offset += cnt;
					cnt = offset;
				}
				for (auto [loc, n] : std::views::reverse(ch_buf)) {
					ch.dat[--cnts_buf[loc]] = n;
				}

				if constexpr (with_planarity) {
					// Set up the reverse mapping for ourselves
					int nxt_ne = ne_en;
					for (int i = ch_en - 1; i >= ch_st; i--) {
						int item = ch.dat[i];
						assert(item >= 1);
						if (item < 1 + NV) continue;
						nxt_ne--;
						rot_edge_ne[item - (1 + NV)] = nxt_ne;
					}
					assert(nxt_ne == ne_st + 1);
					rot_edge_ne[2 * NE] = ne_st;
				}
				auto map_rot_edge = [&] [[gnu::always_inline]] (int ve) -> std::array<int, 4> {
					if constexpr (!with_planarity) return {-1, -1, -1, -1};
					if (!planar) return {-1, -1, -1, -1};
					std::array<int, 4> res{};
					for (int z = 0; z < 4; z++) {
						int o = quarter_edge_matches[4 * ve + z];
						assert(o != -1);
						res[z] = (rot_edge_ne[o >> 2] << 2) + (o & 2) + !(z & 1);
					}
					return res;
				};

				{
					int off = 2 * ne_st;
					for (int i = 2 * nv_st + 1; i <= 2 * nv_en; i++) {
						off += std::exchange(node_adj.bounds[i], off);
					}
					assert(off == 2 * ne_en);
				}

				// Fill in node_edges and node_adj.
				// Reverse order to get the adj in bracket ordering.
				{
					// Handle cap as special: it's first in the node_edges, which means it's in the wrong place for the left endpoint.
					node_adj.bounds[2 * nv_st + 2]++;

					int nxt_ne = ne_en;
					for (int i = ch_en - 1; i >= ch_st; i--) {
						int item = ch.dat[i];
						assert(item >= 1);
						if (item < 1 + NV) continue;
						nxt_ne--;
						auto [v0, v1] = item_vs[item];
						// TODO: Reuse this from the ch pass?
						std::array<int, 2> nvs = {vert_pos_buf[v0], vert_pos_buf[v1]};
						set_ne(nxt_ne, nvs, {
							node_adj.bounds[2 * nvs[0] + 2]++,
							node_adj.bounds[2 * nvs[1] + 1]++,
						}, map_rot_edge(item - (1 + NV)));
					}
					assert(nxt_ne == ne_st + 1);

					// Insert the cap / bump its bound
					set_ne(ne_st, {nv_st, nv_en - 1}, {2 * ne_st, 2 * ne_en - 1}, map_rot_edge(2 * NE));
					node_adj.bounds[2 * nv_en - 1]++;
				}
			} else assert(false);

			int cur_nv = nv_st + (item_vs[cur_item][0] != -1);
			int cur_ne = ne_st + has_cap;
			stk.push_back({cur_idx, ch_st, ch_en, cur_nv, cur_ne});
		};

		auto start_child = [&] [[gnu::always_inline]] () -> int {
			auto& [cur_idx, ch_idx, ch_en, cur_nv, cur_ne] = stk.back();
			assert(ch_idx < ch_en);
			int nxt_item = ch.dat[ch_idx];
			int nxt_idx = nxt_unassigned_idx;
			ch.dat[ch_idx] = nxt_idx;
			par[nxt_idx] = cur_idx;
			int nxt_ne = node_nes.bounds[nxt_idx];
			if (nxt_item < 1 + NV) {
				vert_par_nv[nxt_idx] = cur_nv++;
			} else if (types[cur_idx] != node_type::F && types[cur_idx] != node_type::V) {
				node_edges[cur_ne].twin_ne = nxt_ne;
				node_edges[nxt_ne].twin_ne = cur_ne;
				cur_ne++;
			}

			ch_idx++;
			return nxt_item;
		};

		auto pop_item = [&] [[gnu::always_inline]] () -> void {
			auto [cur_idx, ch_idx, ch_en, cur_nv, cur_ne] = stk.back(); stk.pop_back();
			assert(ch_idx == ch_en);
			subtree_end[cur_idx] = nxt_unassigned_idx;
		};

		par[nxt_unassigned_idx] = -1;
		push_item(ROOT_ITEM);
		while (true) {
			if (stk.back().ch_idx == stk.back().ch_end) {
				pop_item();
				if (stk.empty()) break;
			} else {
				push_item(start_child());
			}
		}

		assert(nxt_unassigned_idx == tot_items);
		assert(ch.bounds.back() == int(ch.dat.size()));
		assert(node_nvs.bounds.back() == int(node_verts.size()));
		assert(node_nes.bounds.back() == int(node_edges.size()));
		assert(node_adj.bounds.back() == int(node_adj.dat.size()));

		// Rewrite node_vertices to the correct index
		for (auto& v : node_verts) {
			v.vert = vert_index[v.vert];
		}

		spqr_tree res{
			std::move(vert_index),
			std::move(edge_index),
			std::move(edge_flipped),
			std::move(par),
			std::move(subtree_end),
			std::move(types),
			std::move(orig_id),
			std::move(ch),
			std::move(node_verts),
			std::move(node_nvs),
			std::move(vert_par_nv),
			std::move(node_edges),
			std::move(node_nes),
			std::move(node_adj),
		};
		if constexpr (with_planarity) {
			return planar_spqr_tree{std::move(res), std::move(node_planar), {std::move(ne_rot_adj)}};
		} else {
			return res;
		}
	}
}

inline std::optional<planar_embedding> planar_embed(const planar_spqr_tree& tree) {
	using node_type = planar_spqr_tree::node_type;

	if (!std::ranges::all_of(tree.node_planar, std::identity{})) {
		return std::nullopt;
	}

	int NE = int(tree.edge_index.size());
	fixed_vector<int> rot_adj(4 * NE, -1);
	auto link = [&] [[gnu::always_inline]] (int a, int b) -> void {
		assert(a != -1 && b != -1);
		assert(rot_adj[a] == -1 && rot_adj[b] == -1);
		assert((a & 1) != (b & 1));
		rot_adj[a] = b;
		rot_adj[b] = a;
	};

	fixed_vector<std::array<std::array<int, 2>, 2>> outer_e(tree.size(), {{{-1, -1}, {-1, -1}}});
	for (int i = tree.size() - 1; i >= 0; i--) {
		auto type = tree.types[i];
		if (type == node_type::F) {
			for (int j : tree.ch[i]) {
				assert(tree.types[j] == node_type::V);
				auto [a, b] = outer_e[j][0];
				if (a != -1) {
					link(a, b);
				}
			}
		} else if (type == node_type::V) {
			std::array<int, 2> qes{-1, -1};
			for (int j : tree.ch[i]) {
				assert(tree.types[j] == node_type::Q);
				if (qes[0] == -1) {
					qes = outer_e[j][0];
				} else {
					link(qes[1], outer_e[j][0][0]);
					qes[1] = outer_e[j][0][1];
				}
			}
			outer_e[i][0] = qes;
		} else if (type == node_type::Q) {
			int e = tree.orig_id[i];
			bool flip = tree.edge_flipped[e];
			std::array<std::array<int, 2>, 2> qes = {{{4 * e + 2 * flip + 0, 4 * e + 2 * flip + 1}, {4 * e + 2 * !flip + 0, 4 * e + 2 * !flip + 1}}};
			if (tree.ch[i].empty()) {
				// Just return ourselves
				outer_e[i] = qes;
			} else {
				int j = tree.ch[i][0];
				if (tree.types[j] == node_type::O) {
					link(qes[0][1], qes[1][0]);
					outer_e[i][0] = {qes[0][0], qes[1][1]};
				} else {
					if (tree.types[j] != node_type::I) {
						link(qes[0][1], outer_e[j][0][0]);
						qes[0][1] = outer_e[j][0][1];
						link(qes[1][0], outer_e[j][1][1]);
						qes[1][0] = outer_e[j][1][0];
					}
					{
						int k = tree.ch[i][1];
						if (outer_e[k][0][0] != -1) {
							link(qes[1][1], outer_e[k][0][0]);
							link(qes[1][0], outer_e[k][0][1]);
						} else {
							link(qes[1][1], qes[1][0]);
						}
					}
					outer_e[i][0] = qes[0];
				}
			}
		} else if (type == node_type::O || type == node_type::I) {
			// Do nothing, the Q node handles it
		} else if (type == node_type::P || type == node_type::S || type == node_type::R) {
			// Just merge things according to the ne_embedding
			for (int ta = 4 * tree.node_nes.bounds[i]; ta < 4 * tree.node_nes.bounds[i+1]; ta++) {
				int tb = tree.ne_embedding.rot_adj[ta];
				if (tb < ta) continue;

				auto tree_qe_to_qe = [&] [[gnu::always_inline]] (int t) -> int {
					return outer_e[tree.node_edges[tree.node_edges[t>>2].twin_ne].node][(t >> 1) & 1][t & 1];
				};
				int qb = tree_qe_to_qe(tb);
				if (ta < 4 * (tree.node_nes.bounds[i] + 1)) {
					outer_e[i][(ta >> 1) & 1][!(ta & 1)] = qb;
				} else {
					int qa = tree_qe_to_qe(ta);
					if ((ta & 3) == 2 && (tb & 3) == 1) {
						// We're the transition between left and right of a vertex, splice it in here.
						int v = tree.node_verts[tree.node_edges[ta >> 2].nvs[1]].vert;
						if (outer_e[v][0][0] != -1) {
							link(qa, outer_e[v][0][1]);
							link(qb, outer_e[v][0][0]);
						} else {
							link(qa, qb);
						}
					} else {
						link(qa, qb);
					}
				}
			}
		} else assert(false);
	}

	return planar_embedding{std::move(rot_adj)};
}

inline std::optional<planar_embedding> planar_embed(
	int NV,
	std::span<const std::array<int, 2>> edges,
	std::span<const int> vert_order,
	std::span<const int> edge_order
) {
	// std::min is by reference, which breaks some optimizations
	auto setmin = [](auto& a, auto b) { if (b < a) a = b; };

	int NE = int(edges.size());
	auto [roots, outedges] = lowval_storted_skeleton_t::build(NV, edges, vert_order, edge_order);

	// Phase 2: do the big ear-decomposition-like walk

	// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
	// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

	// Quarter edges for planar embedding building.
	// Each edge has 4 entries by 4 * edge_id + 2 * source_vert + is_cw (is_cw is arbitrary)
	fixed_vector<int> quarter_edge_matches(4 * NE, -1);
	auto link_quarter_edges = [&] [[gnu::always_inline]] (int a, int b) -> void {
		quarter_edge_matches[a] = b;
		quarter_edge_matches[b] = a;
	};

	{
		int nxt_edge_idx = 0; // Counts backedges only
		fixed_vector<int> postorder_edges(NE);
		auto postorder_edges_end = postorder_edges.begin();
		fixed_vector<bool> postorder_flip(NE+1, false);

		fixed_vector<int> first_occurrence(NV); // First backedge to this depth

		fixed_vector<int> edge_top_depths(NE, -1);

		struct planarity_side_t {
			// For each side, store pointers to the "linked lists" of the edges inside.
			// v[0] is the outer / longer edges and v[1] is the inner / shorter edges, matching the outside-in sort order.
			// The convention is that sides[0].tops[0].depth == top_depth, i.e. at least one minimal return lives on side 0

			// bot_ends are the outer/innermost exposed pieces of the walk down the ear in the tree (they're connected to the bottommost/topmost vertices of the tree path)
			std::array<int, 2> bot_ends{-1, -1};
			struct top_t {
				int end = -1;
				int depth = -1;
			};
			// top_ends are the outer/innermost exposed backedges
			// depths should be increasing going inwards
			std::array<top_t, 2> tops{top_t{-1, -1}, top_t{-1, -1}};
		};
		struct tstack_nonplanarity_t {
			// TODO: What's the nonplanarity certificate look like?
		};
		bounded_stack<planarity_side_t> pstack(with_capacity, NE);
		auto merge_planarity_side = [&] [[gnu::always_inline]] (planarity_side_t& as, const planarity_side_t& bs) -> void {
			assert(as.bot_ends[0] != -1);
			assert(bs.bot_ends[0] != -1);
			link_quarter_edges(as.bot_ends[1], bs.bot_ends[0]);
			as.bot_ends[1] = bs.bot_ends[1];

			if (bs.tops[0].end == -1) {
				// Do nothing
			} else if (as.tops[0].end == -1) {
				as.tops = bs.tops;
			} else {
				// Caller must check that we're planar
				assert(as.tops[1].depth <= bs.tops[0].depth);
				link_quarter_edges(as.tops[1].end, bs.tops[0].end);
				as.tops[1] = bs.tops[1];
			}
		};
		auto make_edge_planarity = [&] [[gnu::always_inline]] (int e_side, int top_depth) -> planarity_side_t {
			edge_top_depths[e_side >> 1] = top_depth;
			return {
				{2 * e_side + 0, 2 * e_side + 1},
				{{{2 * (e_side ^ 1) + 1, top_depth}, {2 * (e_side ^ 1) + 0, top_depth}}},
			};
		};
		struct tstack_t {
			int top_depth;
			int first_idx;
			int pstack_sz;
		};
		bounded_stack<tstack_t> tstack(with_capacity, NV + NE);
		auto cur_tstack = [&] [[gnu::always_inline]] () -> tstack_t& { return tstack.end()[-1]; };
		auto nxt_tstack = [&] [[gnu::always_inline]] () -> tstack_t& { return tstack.end()[-2]; };

		auto push_tstack = [&] [[gnu::always_inline]] (int top_depth, int pstack_sz) -> void {
			tstack.emplace_back(top_depth, nxt_edge_idx, pstack_sz);
		};
		auto push_vert_tstack = [&] [[gnu::always_inline]] (int top_depth) -> void {
			push_tstack(top_depth, 0);
		};
		auto push_edge_tstack = [&] [[gnu::always_inline]] (int top_depth, int e_side) -> int {
			pstack.emplace_back(make_edge_planarity(e_side, top_depth));
			push_tstack(top_depth, 1);
			*postorder_edges_end++ = (e_side >> 1);
			return nxt_edge_idx++;
		};
		auto flip_tstack_planarity = [&] [[gnu::always_inline]] (tstack_t* t) -> void {
			postorder_flip[t->first_idx] ^= 1;
			postorder_flip[(t+1==tstack.end()) ? nxt_edge_idx : (t+1)->first_idx] ^= 1;
		};

		struct dfs_stack_t {
			bool has_vert_tstack;
			int ch_idx;
			int ch_end;
			int orig_tstack;
		};
		bounded_stack<dfs_stack_t> stk(with_capacity, NV);
		for (auto rt : roots) {
			auto push_vert = [&] [[gnu::always_inline]] (int cur) -> void {
				int cur_depth = int(stk.size());

				int lo = outedges.bounds[cur];
				int hi = outedges.bounds[cur+1];
				bool has_vert_tstack;
				{
					// Find the first same-BCC edge, and check it's type 2 (has lowval2), if so it's the ear tstack and we defer pushing ourselves.
					int first_edge = lo;
					while (first_edge < hi && outedges.dat[first_edge].key.is_new_block()) first_edge++;
					if (first_edge < hi && outedges.dat[first_edge].key.is_type_2()) {
						// Move first_edge to the beginning
						auto e = outedges.dat[first_edge];
						std::move_backward(outedges.dat.begin() + lo, outedges.dat.begin() + first_edge, outedges.dat.begin() + first_edge + 1);
						outedges.dat[lo] = e;
						has_vert_tstack = false;
					} else {
						push_vert_tstack(cur_depth);
						has_vert_tstack = true;
					}
				}
				stk.push_back({has_vert_tstack, lo, hi, -1});
			};
			// return true means jump to start_edge, return false means jump to finish_edge
			auto start_edge = [&] [[gnu::always_inline]] () -> std::optional<int> {
				int cur_depth = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);
				auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
				auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

				if (lowval >= cur_depth || is_type_1) assert(s.has_vert_tstack);

				s.orig_tstack = int(tstack.size());
				if (is_tree) {
					first_occurrence[cur_depth] = NE;
					return nxt;
				} else {
					return std::nullopt;
				}
			};
			auto finish_edge = [&] [[gnu::always_inline]] [[nodiscard]] () -> std::optional<tstack_nonplanarity_t> {
				int cur_depth = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);

				auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
				s.ch_idx++;

				auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

				const auto orig_tstack_end = tstack.begin() + s.orig_tstack;

				auto merge_v_s = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 1);
					assert(nxt_tstack().pstack_sz <= 1);
					if (nxt_tstack().pstack_sz) {
						// TODO: Can inline further
						merge_planarity_side(pstack.end()[-2], pstack.end()[-1]);
						pstack.pop_back();
					}
					nxt_tstack().top_depth = cur_tstack().top_depth;
					nxt_tstack().pstack_sz = 1;
					tstack.pop_back();
				};

				auto merge_v_r = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 2);
					assert(nxt_tstack().pstack_sz <= 1);
					if (nxt_tstack().pstack_sz) {
						// TODO: Can inline further
						merge_planarity_side(pstack.end()[-3], pstack.end()[-2]);
						pstack.end()[-2] = pstack.end()[-1];
						pstack.pop_back();
					}
					nxt_tstack().top_depth = cur_tstack().top_depth;
					nxt_tstack().pstack_sz = 2;
					tstack.pop_back();
				};

				auto merge_s_s = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 1);
					assert(nxt_tstack().pstack_sz == 1);
					auto& a = pstack.end()[-2];
					auto& b = pstack.end()[-1];
					link_quarter_edges(a.tops[0].end, b.bot_ends[0]);
					link_quarter_edges(a.tops[1].end, b.bot_ends[1]);
					a.tops = b.tops;
					pstack.pop_back();
					nxt_tstack().top_depth = cur_depth;
					tstack.pop_back();
				};

				auto merge_r_s_and_fold = [&] [[gnu::always_inline]] () -> void {
					assert(nxt_tstack().pstack_sz == 2);
					// First, just merge cur_tstack() as a backedge
					merge_planarity_side(pstack.end()[-3], pstack.end()[-1]);
					pstack.pop_back();
					// Fixup: flatten the R-node into a single backedge
					auto& s0 = pstack.end()[-2];
					auto& s1 = pstack.end()[-1];
					link_quarter_edges(s0.bot_ends[1], s1.bot_ends[1]);
					s0.bot_ends[1] = s1.bot_ends[0];
					if (s1.tops[1].end != -1) {
						link_quarter_edges(s0.tops[1].end, s1.tops[1].end);
						s0.tops[1] = s1.tops[0];
					}
					pstack.pop_back();
					nxt_tstack().pstack_sz = 1;
					tstack.pop_back();
				};


				auto merge_p_s = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 1);
					assert(nxt_tstack().pstack_sz == 1);
					// TODO: Can inline further
					merge_planarity_side(pstack.end()[-2], pstack.end()[-1]);
					pstack.pop_back();
					tstack.pop_back();
				};

				auto merge_s_r = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 2);
					assert(nxt_tstack().pstack_sz == 1);
					{
						auto& c = pstack.end()[-3];
						auto& s0 = pstack.end()[-2];
						auto& s1 = pstack.end()[-1];
						link_quarter_edges(c.tops[0].end, s0.bot_ends[0]);
						link_quarter_edges(c.tops[1].end, s1.bot_ends[0]);
						s0.bot_ends[0] = c.bot_ends[0];
						s1.bot_ends[0] = c.bot_ends[1];
						c = s0;
						s0 = s1;
						pstack.pop_back();
					}
					nxt_tstack().top_depth = cur_tstack().top_depth;
					nxt_tstack().pstack_sz = 2;
					tstack.pop_back();
				};

				auto merge_r_r = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 2);
					assert(nxt_tstack().pstack_sz == 2);
					merge_planarity_side(pstack.end()[-3], pstack.end()[-1]);
					merge_planarity_side(pstack.end()[-4], pstack.end()[-2]);
					pstack.pop_back();
					pstack.pop_back();
					tstack.pop_back();
				};

				auto merge_p_r_side_0 = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 2);
					assert(nxt_tstack().pstack_sz == 1);
					merge_planarity_side(pstack.end()[-3], pstack.end()[-2]);
					pstack.end()[-2] = pstack.end()[-1];
					pstack.pop_back();
					nxt_tstack().pstack_sz = 2;
					tstack.pop_back();
				};

				auto merge_p_r_side_1 = [&] [[gnu::always_inline]] () -> void {
					assert(cur_tstack().pstack_sz == 2);
					assert(nxt_tstack().pstack_sz == 1);
					merge_planarity_side(pstack.end()[-3], pstack.end()[-1]);
					std::swap(pstack.end()[-3], pstack.end()[-2]);
					pstack.pop_back();
					nxt_tstack().pstack_sz = 2;
					tstack.pop_back();
				};

				if (lowval >= cur_depth) {
					if (is_tree) {
						push_edge_tstack(cur_depth, e_side ^ 1);
						if (lowval == cur_depth) {
							// Merge the backedge
							merge_p_s();
						}
						// Merge the vertex
						merge_v_s();
						{
							// Join the bottom together
							auto& p = pstack.back();
							link_quarter_edges(p.bot_ends[0], p.bot_ends[1]);
							p.bot_ends[0] = p.tops[1].end;
							p.bot_ends[1] = p.tops[0].end;
							p.tops = {};
						}
					} else {
						push_edge_tstack(lowval, e_side);
						assert(cur_tstack().pstack_sz == 1);
						{
							// Join the loop together
							auto& p = pstack.back();
							link_quarter_edges(p.bot_ends[1], p.tops[1].end);
							p.bot_ends[1] = p.tops[0].end;
							p.tops = {};
						}
					}
					assert(s.has_vert_tstack);
					// Merge into the vertex tstack
					merge_v_s();
					return std::nullopt;
				}

				if (is_tree) {
					push_edge_tstack(cur_depth, e_side ^ 1);
					while (nxt_tstack().top_depth >= cur_depth) {
						if (nxt_tstack().top_depth > cur_depth) {
							if (tstack.end()[-3].top_depth < cur_depth) {
								break;
							}
							// Merge the vertex in
							merge_v_s();
							if (nxt_tstack().top_depth > cur_depth) {
								// S-type merge
								merge_s_s();
							} else {
								assert(nxt_tstack().top_depth == cur_depth);
								// R-type merge
								merge_r_s_and_fold();
							}
						} else {
							// P-type merge
							merge_p_s();
						}
					}

					if (cur_tstack().first_idx > first_occurrence[cur_depth]) {
						auto source = tstack.end() - 2;
						while (source->first_idx > first_occurrence[cur_depth]) --source;

						// We won't care about bot_ends[1] here at all
						// TODO: This shouldn't live on the pstack
						pstack.push_back({});
						pstack.end()[-1].bot_ends[0] = std::exchange(pstack.end()[-2].bot_ends[1], -1);
						cur_tstack().pstack_sz = 2;

						// last_top == cur_tstack().top_depth
						int last_top = cur_depth;
						while (true) {
							if (nxt_tstack().top_depth > cur_depth) {
								// TODO: If we have separate atoms, coalesce them now
								// Vertex, check the edge instead
								merge_v_r();
								if (nxt_tstack().top_depth > cur_depth) {
									// Tree edge, always fine, just extend
									assert(tstack.end() > source + 2);
									// last_top < cur_depth, since otherwise we would've merged above
									assert(last_top < cur_depth);
									merge_s_r();
								} else {
									// Chunk entry
									// TODO: If we have separate atoms, handle this correctly
									assert(nxt_tstack().pstack_sz == 2);
									// -4 is side 0, -3 is side 1
									if (tstack.end() == source + 2) {
										// TODO: separate atom handling
										int t0 = pstack.end()[-4].tops[1].depth;
										int t1 = pstack.end()[-3].tops[1].depth;
										assert(t0 == cur_depth || t1 == cur_depth);
										assert(t0 != -1);
										if (std::min(t0, t1) > last_top) {
											assert(last_top < cur_depth);
											return tstack_nonplanarity_t{};
										}
										if (t1 != cur_depth) {
											flip_tstack_planarity(tstack.end() - 1);
											std::swap(pstack.end()[-2], pstack.end()[-1]);
											merge_r_r();
											if (last_top < cur_tstack().top_depth) {
												std::swap(pstack.end()[-2], pstack.end()[-1]);
												flip_tstack_planarity(tstack.end() - 1);
												cur_tstack().top_depth = last_top;
											}
										} else {
											merge_r_r();
										}
										break;
									}
									if (nxt_tstack().top_depth == cur_depth) {
										assert(last_top < cur_depth);
										// Never have any atoms here
										assert(nxt_tstack().pstack_sz == 2);
										if (pstack.end()[-3].tops[0].end != -1) {
											// Double-sided to cur_depth, conflicts with cur_tstack()
											return tstack_nonplanarity_t{};
										}
										// Flip: we will put cur_depth on side 1 until the bottom
										std::swap(pstack.end()[-4], pstack.end()[-3]);
										flip_tstack_planarity(tstack.end() - 2);
										nxt_tstack().top_depth = last_top;
									} else {
										// TODO: With atoms, we should check pstack_idx+1
										if (pstack.end()[-3].tops[0].end != -1 && pstack.end()[-3].tops[0].depth != cur_depth) {
											// Non-empty on both sides, conflicts with source
											return tstack_nonplanarity_t{};
										}
										if (pstack.end()[-4].tops[1].depth > last_top) {
											// 3 nonlaminar edges with nxt_tstack(), cur_tstack(), source
											return tstack_nonplanarity_t{};
										}
										last_top = nxt_tstack().top_depth;
									}
									merge_r_r();
								}
							} else {
								// Single atom
								assert(nxt_tstack().pstack_sz == 1);
								if (nxt_tstack().top_depth == cur_depth || tstack.end() == source + 2) {
									// We will put cur_depth on side 1 until the bottom
									flip_tstack_planarity(tstack.end() - 2);
									merge_p_r_side_1();
									if (tstack.end() == source + 1) {
										if (cur_tstack().top_depth <= last_top) {
											// Flip it back
											std::swap(pstack.end()[-2], pstack.end()[-1]);
											flip_tstack_planarity(tstack.end() - 1);
										} else {
											cur_tstack().top_depth = last_top;
										}
										break;
									}
									cur_tstack().top_depth = last_top;
								} else {
									if (pstack.end()[-3].tops[1].depth > last_top) {
										// 3 nonlaminar edges with nxt_tstack(), cur_tstack(), source
										return tstack_nonplanarity_t{};
									}
									last_top = nxt_tstack().top_depth;
									// Merge into side 0
									merge_p_r_side_0();
								}
							}
						}

						assert(tstack.end() == source + 1);
						assert(cur_tstack().top_depth < cur_depth);

						// Link inner ones
						link_quarter_edges(pstack.end()[-2].tops[1].end, pstack.end()[-1].tops[1].end);

						// Prune cur_depth things from tops of each end
						// TODO: Handle atoms correctly
						for (auto& side : std::span(pstack.end() - 2, pstack.end())) {
							assert(side.tops[1].depth == cur_depth);
							auto t = side.tops[1].end;
							while (true) {
								int nt = quarter_edge_matches[t ^ 1];
								if (nt == -1) {
									side.bot_ends[1] = t^1;
									side.tops = {};
									break;
								}
								int nd = edge_top_depths[nt >> 2];
								assert(nd <= cur_depth);
								if (nd < cur_depth) {
									quarter_edge_matches[nt] = -1;
									quarter_edge_matches[t^1] = -1;
									side.bot_ends[1] = t^1;
									side.tops[1] = {nt, nd};
									break;
								}
								t = nt;
							}
						}
					}

					if (is_type_1) assert(s.has_vert_tstack);
					if (s.has_vert_tstack) {
						// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
						assert(tstack.end() >= orig_tstack_end + 3);

						if (cur_tstack().top_depth == cur_depth) {
							// We're currently a backedge, expand to a 2-sided chunk
							assert(cur_tstack().pstack_sz == 1);
							pstack.push_back({});
							auto& s0 = pstack.end()[-2];
							auto& s1 = pstack.end()[-1];
							s1.bot_ends[0] = s0.bot_ends[1];
							s0.bot_ends[1] = s0.tops[0].end;
							s1.bot_ends[1] = s0.tops[1].end;
							s0.tops = {};
							cur_tstack().pstack_sz = 2;
						}

						if (!is_type_1) {
							// The lowval side should be side 1, everything else goes on side 0.
							// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
							// but is guaranteed to have *something* > lowval by non-type-1-ness
							int last_top = cur_depth;
							for (auto t = tstack.end(); t >= orig_tstack_end + 2; --t) {
								if (t == tstack.end() || t->top_depth >= cur_depth) {
									// We're a vertex
									if (t != tstack.end()) {
										assert(t == tstack.end() - 2);
										merge_v_r();
									}
									--t;
									if (t->top_depth >= cur_depth) {
										// S-type merge, nothing can break
										if (t != tstack.end() - 1) {
											assert(t == tstack.end() - 2);
											merge_s_r();
										}
									} else {
										// R-type chunk
										assert(t->pstack_sz == 2);
										auto& s0 = pstack.end()[-2 * (t == tstack.end() - 2) - 2];
										auto& s1 = pstack.end()[-2 * (t == tstack.end() - 2) - 1];
										if (t->top_depth == lowval) {
											if (t > orig_tstack_end + 2 || s0.tops[1].depth == lowval) {
												flip_tstack_planarity(t);
												std::swap(s0, s1);
											}
										}
										int next_top = s0.tops[1].depth;
										if (t == orig_tstack_end + 2) assert(next_top != -1 && next_top > lowval);
										if (s1.tops[1].end != -1 && s1.tops[1].depth != lowval) {
											return tstack_nonplanarity_t{};
										}
										if (next_top != -1) {
											if (next_top > last_top) return tstack_nonplanarity_t{};
											last_top = s0.tops[0].depth;
										}

										if (t != tstack.end() - 1) {
											assert(t == tstack.end() - 2);
											merge_r_r();
										}
									}
								} else {
									assert(t == tstack.end() - 2);
									// single atom
									assert(t->pstack_sz == 1);
									if (t > orig_tstack_end + 2 && t->top_depth == lowval) {
										if (pstack.end()[-3].tops[1].depth > lowval) {
											return tstack_nonplanarity_t{};
										}
										// Merge into side 1
										flip_tstack_planarity(t);
										merge_p_r_side_1();
									} else {
										int next_top = pstack.end()[-3].tops[1].depth;
										assert(next_top != -1);
										// Either i == orig_tstack + 2, or we flipped already
										assert(next_top > lowval);
										if (next_top > last_top) return tstack_nonplanarity_t{};

										// Just merge into side 0
										last_top = t->top_depth;
										merge_p_r_side_0();
									}
								}
							}
						}

						assert(tstack.end() == orig_tstack_end + 3);

						// Merge with the backedge
						merge_p_r_side_0();
						// Merge with the vertex
						merge_v_r();

						assert(cur_tstack().top_depth == lowval);

						{
							// Join the 2 sides to 1 big backedge
							assert(cur_tstack().pstack_sz == 2);
							auto& s0 = pstack.end()[-2];
							auto& s1 = pstack.end()[-1];
							assert(s1.bot_ends[1] != -1);
							link_quarter_edges(s0.bot_ends[0], s1.bot_ends[0]);
							s0.bot_ends[0] = s1.bot_ends[1];
							if (s1.tops[0].end != -1) {
								assert(s1.tops[1].depth == lowval);
								assert(s1.tops[0].depth == lowval);
								link_quarter_edges(s0.tops[0].end, s1.tops[0].end);
								s0.tops[0].end = s1.tops[1].end;
								// Already true since the backedge was on side 0
								assert(s0.tops[0].depth == lowval);
							}
							pstack.pop_back();
							cur_tstack().pstack_sz = 1;
						}
						assert(pstack.end()[-1].bot_ends[1] != -1);
						assert(pstack.end()[-1].bot_ends[0] != -1);
					}
				} else {
					assert(is_type_1);
					int idx = push_edge_tstack(lowval, e_side);
					setmin(first_occurrence[lowval], idx);
				}

				if (is_type_1 && nxt_tstack().top_depth == lowval) {
					assert(s.has_vert_tstack);
					merge_p_s();
				}

				if (!s.has_vert_tstack) {
					assert(!is_type_1);
					// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
					push_vert_tstack(cur_depth);
					s.has_vert_tstack = true;
				}
				return std::nullopt;
			};
			auto pop_vert = [&] [[gnu::always_inline]] () -> void {
				auto& s = stk.back();
				assert(s.ch_idx == s.ch_end);
				assert(s.has_vert_tstack);
				stk.pop_back();
			};

			push_vert(rt);
			while (true) {
				if (stk.back().ch_idx == stk.back().ch_end) {
					pop_vert();
					if (stk.empty()) break;
					if (auto res = finish_edge(); res) return std::nullopt;
				} else if (std::optional<int> nxt = start_edge(); nxt) {
					push_vert(*nxt);
				} else {
					if (auto res = finish_edge(); res) assert(false);
				}
			}
			assert(int(tstack.size()) == 1);
			assert(int(pstack.size()) <= 1);
			if (!pstack.empty()) {
				auto& s0 = pstack.back();
				// Fold the root
				int a = s0.bot_ends[0];
				int b = s0.bot_ends[1];
				quarter_edge_matches[a] = b;
				quarter_edge_matches[b] = a;
				pstack.pop_back();
			}
			tstack.pop_back();
		}
		assert(nxt_edge_idx == NE);
		{
			fixed_vector<bool> edge_flip(NE, false);
			{
				bool planarity_flip = false;
				for (int e = 0; e < NE; e++) {
					planarity_flip ^= postorder_flip[e];
					edge_flip[postorder_edges[e]] = planarity_flip;
				}
				planarity_flip ^= postorder_flip[NE];
				assert(!planarity_flip);
			}
			for (int e = 0; e < NE; e++) {
				if (edge_flip[e]) {
					std::swap(quarter_edge_matches[4*e + 0], quarter_edge_matches[4*e + 1]);
					std::swap(quarter_edge_matches[4*e + 2], quarter_edge_matches[4*e + 3]);
				}
				for (int z = 0; z < 4; z++) {
					quarter_edge_matches[4*e + z] = (quarter_edge_matches[4*e+z] >> 1 << 1) | !(z & 1);
				}
			}
		}
	}
	return planar_embedding{std::move(quarter_edge_matches)};
}

inline bool can_planar_embed(
	int NV,
	std::span<const std::array<int, 2>> edges,
	std::span<const int> vert_order,
	std::span<const int> edge_order
) {
	// std::min is by reference, which breaks some optimizations
	auto setmin = [](auto& a, auto b) { if (b < a) a = b; };

	int NE = int(edges.size());

	auto [roots, outedges] = lowval_storted_skeleton_t::build(NV, edges, vert_order, edge_order);

	// Phase 2: do the big ear-decomposition-like walk

	// We're going to build a tree of all SPQR *nodes* + all original *vertices* (collectively *items*).
	// Vertices will hang off the first SPQR node containing them, and blocks will be rooted at a topmost Q node for the top edge.

	{
		int nxt_edge_idx = 0; // Counts backedges only

		fixed_vector<int> first_occurrence(NV); // First backedge to this depth

		struct tstack_planarity_side_t {
			// For each side, store pointers to the "linked lists" of the edges inside.
			// v[0] is the outer / longer edges and v[1] is the inner / shorter edges, matching the outside-in sort order.

			struct top_t {
				int end = -1;
				int depth = -1;
			};
			// top_ends are the outer/innermost exposed backedges
			// depths should be increasing going inwards
			std::array<top_t, 2> tops{top_t{-1, -1}, top_t{-1, -1}};
		};
		struct tstack_planarity_t {
			// The convention is that sides[0].tops[0].depth == top_depth, i.e. at least one minimal return lives on side 0
			std::array<tstack_planarity_side_t, 2> sides;
		};
		fixed_vector<tstack_planarity_side_t::top_t> prev_edge(NE, {-1, -1});
		auto merge_planarity_side = [&] [[gnu::always_inline]] (tstack_planarity_side_t& as, const tstack_planarity_side_t& bs) -> void {
			// If there's no bottom edges, then we must be an isolated vertex, so we can end early.
			// Caller must check that we're planar
			assert(as.tops[1].depth <= bs.tops[0].depth);
			prev_edge[bs.tops[0].end] = as.tops[1];
			as.tops[1] = bs.tops[1];
		};
		auto make_edge_planarity = [&] [[gnu::always_inline]] (int e, int top_depth, bool is_tree) -> tstack_planarity_t {
			tstack_planarity_t p;
			if (!is_tree) {
				p.sides[0].tops = {{{e, top_depth}, {e, top_depth}}};
			}
			return p;
		};
		struct tstack_t {
			int top_depth = -1;
			int first_idx = -1;
			tstack_planarity_t planarity;
		};
		bounded_stack<tstack_t> tstack(with_capacity, NV + NE);
		auto cur_tstack = [&] [[gnu::always_inline]] () -> tstack_t& { return tstack.end()[-1]; };
		auto nxt_tstack = [&] [[gnu::always_inline]] () -> tstack_t& { return tstack.end()[-2]; };

		auto push_tstack = [&] [[gnu::always_inline]] (int top_depth, tstack_planarity_t planarity) -> void {
			tstack.emplace_back(top_depth, nxt_edge_idx, planarity);
		};
		auto push_vert_tstack = [&] [[gnu::always_inline]] (int top_depth) -> void {
			push_tstack(top_depth, {});
		};
		auto push_edge_tstack = [&] [[gnu::always_inline]] (int top_depth, int e, bool is_tree) -> int {
			push_tstack(top_depth, make_edge_planarity(e, top_depth, is_tree));
			return nxt_edge_idx++;
		};

		struct dfs_stack_t {
			bool has_vert_tstack;
			int ch_idx;
			int ch_end;
			int orig_tstack;
		};
		bounded_stack<dfs_stack_t> stk(with_capacity, NV);
		for (auto rt : roots) {
			auto push_vert = [&] [[gnu::always_inline]] (int cur) -> void {
				int cur_depth = int(stk.size());

				int lo = outedges.bounds[cur];
				int hi = outedges.bounds[cur+1];
				bool has_vert_tstack;
				{
					// Find the first same-BCC edge, and check it's type 2 (has lowval2), if so it's the ear tstack and we defer pushing ourselves.
					int first_edge = lo;
					while (first_edge < hi && outedges.dat[first_edge].key.is_new_block()) first_edge++;
					if (first_edge < hi && outedges.dat[first_edge].key.is_type_2()) {
						// Move first_edge to the beginning
						auto e = outedges.dat[first_edge];
						std::move_backward(outedges.dat.begin() + lo, outedges.dat.begin() + first_edge, outedges.dat.begin() + first_edge + 1);
						outedges.dat[lo] = e;
						has_vert_tstack = false;
					} else {
						push_vert_tstack(cur_depth);
						has_vert_tstack = true;
					}
				}
				stk.push_back({has_vert_tstack, lo, hi, -1});
			};
			// return true means jump to start_edge, return false means jump to finish_edge
			auto start_edge = [&] [[gnu::always_inline]] () -> std::optional<int> {
				int cur_depth = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);
				auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
				auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

				if (lowval >= cur_depth || is_type_1) assert(s.has_vert_tstack);

				s.orig_tstack = int(tstack.size());
				if (is_tree) {
					first_occurrence[cur_depth] = NE;
					return nxt;
				} else {
					return std::nullopt;
				}
			};
			auto finish_edge = [&] [[gnu::always_inline]] [[nodiscard]] () -> bool {
				int cur_depth = int(stk.size()) - 1;
				auto& s = stk.back();
				assert(s.ch_idx < s.ch_end);

				auto [_, nxt, e_side, key] = outedges.dat[s.ch_idx];
				int e = e_side >> 1;
				s.ch_idx++;

				auto [lowval, is_tree, is_type_1] = key.unpack(cur_depth);

				const int orig_tstack = s.orig_tstack;

				if (lowval >= cur_depth) {
					if (is_tree) {
						if (lowval == cur_depth) {
							// Delete the backedge
							tstack.pop_back();
						}
						// Delete the vertex
						tstack.pop_back();
					} else {
					}
					assert(s.has_vert_tstack);
					assert(int(tstack.size()) == orig_tstack);
					return true;
				}

				if (is_tree) {
					push_edge_tstack(cur_depth, e, true);
					while (nxt_tstack().top_depth >= cur_depth) {
						if (nxt_tstack().top_depth > cur_depth) {
							if (tstack.end()[-3].top_depth < cur_depth) {
								break;
							}
							// Merge the vertex in
							tstack.pop_back();
						}

						tstack.pop_back();
					}
					cur_tstack().planarity = tstack_planarity_t{};

					if (cur_tstack().first_idx > first_occurrence[cur_depth]) {
						int source = int(tstack.size()) - 2;
						while (tstack[source].first_idx > first_occurrence[cur_depth]) --source;

						// last_top == cur_tstack().top_depth
						int last_top = cur_depth;
						tstack_planarity_side_t cur_planarity = {};
						while (int(tstack.size()) > source + 2) {
							if (nxt_tstack().top_depth > cur_depth) {
								// Vertex or tree edge, no conditions
							} else if (nxt_tstack().top_depth == cur_depth) {
								if (nxt_tstack().planarity.sides[1].tops[0].depth != -1) {
									// Double-sided to cur_depth, conflicts with cur_tstack()
									assert(last_top < cur_depth);
									return false;
								}
								// Throw away the inner edge
							} else {
								if (nxt_tstack().planarity.sides[1].tops[0].depth != -1 && nxt_tstack().planarity.sides[1].tops[0].depth != cur_depth) {
									// Non-empty on both sides, conflicts with source
									return false;
								}
								if (nxt_tstack().planarity.sides[0].tops[1].depth > last_top) {
									// Nonlaminar with cur_tstack()
									return false;
								}
								if (last_top < cur_depth) {
									auto nxt_planarity = nxt_tstack().planarity.sides[0];
									prev_edge[cur_planarity.tops[0].end] = nxt_planarity.tops[1];
									cur_planarity.tops[0] = nxt_planarity.tops[0];
								} else {
									cur_planarity = nxt_tstack().planarity.sides[0];
								}
								last_top = nxt_tstack().top_depth;
							}
							tstack.pop_back();
						}

						int t0 = nxt_tstack().planarity.sides[0].tops[1].depth;
						int t1 = nxt_tstack().planarity.sides[1].tops[1].depth;
						assert(t0 == cur_depth || t1 == cur_depth);
						// Handles -1 correctly
						if (std::min(t0, t1) > last_top) {
							assert(last_top < cur_depth);
							return false;
						}
						if (last_top < cur_depth) {
							if (t0 == cur_depth) {
								// We need to flip cur_tstack and nxt_tstack relative to each other.
								// Flip the one with worse top_depth.
								if (t1 != -1) {
									// merge into side 1
									merge_planarity_side(nxt_tstack().planarity.sides[1], cur_planarity);
								} else {
									nxt_tstack().planarity.sides[1] = cur_planarity;
								}
								if (last_top < nxt_tstack().top_depth) {
									nxt_tstack().top_depth = last_top;
									std::swap(nxt_tstack().planarity.sides[0], nxt_tstack().planarity.sides[1]);
								}
							} else {
								assert(t0 < cur_depth);
								assert(t0 <= last_top);
								// merge into side 0
								merge_planarity_side(nxt_tstack().planarity.sides[0], cur_planarity);
							}
						}
						tstack.pop_back();

						// Prune off finished cur-side things
						for (auto& side : cur_tstack().planarity.sides) {
							while (side.tops[1].depth == cur_depth) {
								side.tops[1] = prev_edge[side.tops[1].end];
								if (side.tops[1].end == -1) {
									side.tops[0] = {};
								}
							}
						}
					}

					if (is_type_1) assert(s.has_vert_tstack);
					if (s.has_vert_tstack) {
						// NB: tstack[orig_size] is the vertex and tstack[orig_size+1] is the backedge; maybe we should reverse them?
						assert(int(tstack.size()) >= orig_tstack + 3);

						if (!is_type_1) {
							// The lowval side should be side 1, everything else goes on side 0.
							// The exception is tstack[orig_tstack + 2], which could be == lowval on one/both sides,
							// but is guaranteed to have *something* > lowval by non-type-1-ness
							auto& t = tstack[orig_tstack + 2];
							{
								assert(t.planarity.sides[0].tops[0].depth == t.top_depth);
								assert(t.planarity.sides[0].tops[1].depth != -1);
								if (t.planarity.sides[0].tops[1].depth == lowval) {
									t.planarity.sides[0] = t.planarity.sides[1];
								} else if (t.planarity.sides[1].tops[1].depth != -1 && t.planarity.sides[1].tops[1].depth != lowval) {
									return false;
								}
								assert(t.planarity.sides[0].tops[1].depth > lowval);
							}
							tstack_planarity_side_t cur_planarity = t.planarity.sides[0];
							for (int i = orig_tstack + 3; i < int(tstack.size()); i++) {
								if (tstack[i].top_depth == lowval) {
									std::swap(tstack[i].planarity.sides[0], tstack[i].planarity.sides[1]);
								}
								if (tstack[i].planarity.sides[1].tops[1].depth != -1 && tstack[i].planarity.sides[1].tops[1].depth != lowval) {
									return false;
								}
								int next_top = tstack[i].planarity.sides[0].tops[0].depth;
								if (next_top != -1) {
									if (cur_planarity.tops[1].depth > next_top) return false;
									merge_planarity_side(cur_planarity, tstack[i].planarity.sides[0]);
								}
							}
							merge_planarity_side(tstack[orig_tstack+1].planarity.sides[0], cur_planarity);
							tstack[orig_tstack+1].planarity.sides[1] = tstack_planarity_side_t{};
						} else {
							assert(int(tstack.size()) == orig_tstack + 3);
						}
						tstack[orig_tstack] = tstack[orig_tstack+1];
						tstack.truncate(orig_tstack + 1);
						assert(cur_tstack().top_depth == lowval);
					}
				} else {
					assert(is_type_1);
					int idx = push_edge_tstack(lowval, e, false);
					setmin(first_occurrence[lowval], idx);
				}

				if (is_type_1 && nxt_tstack().top_depth == lowval) {
					assert(s.has_vert_tstack);
					tstack.pop_back();
				}

				if (!s.has_vert_tstack) {
					assert(!is_type_1);
					// Throw cur_vert_node onto the tstack so it'll get interleaved correctly
					push_vert_tstack(cur_depth);
					s.has_vert_tstack = true;
				}
				return true;
			};
			auto pop_vert = [&] [[gnu::always_inline]] () -> void {
				auto& s = stk.back();
				assert(s.ch_idx == s.ch_end);
				assert(s.has_vert_tstack);
				stk.pop_back();
			};

			push_vert(rt);
			while (true) {
				if (stk.back().ch_idx == stk.back().ch_end) {
					pop_vert();
					if (stk.empty()) break;
					if (auto res = finish_edge(); !res) return false;
				} else if (std::optional<int> nxt = start_edge(); nxt) {
					push_vert(*nxt);
				} else {
					if (auto res = finish_edge(); !res) assert(false);
				}
			}
			assert(int(tstack.size()) == 1);
			tstack.pop_back();
		}
	}
	return true;
}

} // namespace wala
