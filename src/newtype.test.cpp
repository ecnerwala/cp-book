#include <catch2/catch_test_macros.hpp>

#include "newtype.hpp"

#include <cstdint>
#include <numeric>
#include <type_traits>

INT_NEWTYPE(node_t);
INT_NEWTYPE(edge_t);
INT_NEWTYPE(big_t, int64_t);

// Distinct tags give distinct, non-interconvertible types
static_assert(!std::is_same_v<node_t, edge_t>);
static_assert(!std::is_convertible_v<node_t, edge_t>);
static_assert(!std::is_convertible_v<edge_t, node_t>);
static_assert(std::is_constructible_v<edge_t, node_t>);
static_assert(std::is_constructible_v<big_t, node_t>);

// Raw ints embed implicitly (exact underlying type at runtime, any integral
// constant at compile time); conversion back out is explicit
static_assert(std::is_convertible_v<int, node_t>);
static_assert(!std::is_convertible_v<node_t, int>);
// (runtime conversion from other integral types, e.g. int64_t -> node_t, goes
// through the consteval ctor and is rejected unless it's a constant expression)

TEST_CASE("Newtype basics", "[newtype]") {
	node_t a = 3;   // implicit from literal (consteval)
	node_t b(7);    // also implicit from runtime ints of the underlying type
	REQUIRE(a != b);
	REQUIRE(a < b);
	REQUIRE(a.v == 3);
	REQUIRE(int(b) == 7);

	// ring arithmetic, closed within the type; literals participate implicitly
	node_t c = a + 4;
	REQUIRE(c == b);
	REQUIRE(b - a == node_t(4));
	REQUIRE(bool(b - a == 4)); // bare literals work outside Catch's expression decomposer
	c += 2;
	REQUIRE(c - b == node_t(2));
	REQUIRE(--c == b + 1);
	REQUIRE(c++ == b + 1);
	REQUIRE(c == b + 2);

	// 0 and 1 are distinguished points; *, /, %, shifts are closed too
	REQUIRE(a * 0 == node_t(0));
	REQUIRE(b * 1 == b);
	REQUIRE(a * b == node_t(21));
	REQUIRE(b / 2 == a);
	REQUIRE(b % a == node_t(1));
	REQUIRE((a << 1) == node_t(6));
	REQUIRE((b >> 1) == a);
	REQUIRE((a & b) == a);
	REQUIRE(-a == node_t(-3));
	// a + node_t2(1) does not compile: cross-type operands need explicit casts

	// explicit casts between newtypes
	edge_t e(a);
	REQUIRE(e == edge_t(3));
	big_t big(b);
	REQUIRE(big.v == int64_t(7));
}

TEST_CASE("Typed vector", "[newtype]") {
	INT_NEWTYPE(idx_t);
	vec<idx_t, int> a(idx_t(4), 10);
	REQUIRE(a.size() == idx_t(4));
	a[idx_t(2)] = 5;
	REQUIRE(a[idx_t(2)] == 5);
	REQUIRE(a.at(idx_t(0)) == 10);

	idx_t k = a.push_back(99);
	REQUIRE(k == idx_t(4));
	REQUIRE(a.back() == 99);

	int sum = 0;
	for (int x : a) sum += x;
	REQUIRE(sum == 10 + 10 + 5 + 10 + 99);

	int sum2 = 0;
	for (idx_t i = 0; i < a.size(); ++i) sum2 += a[i];
	REQUIRE(sum2 == sum);

	int sum3 = 0;
	for (idx_t i : a.keys()) sum3 += a[i];
	REQUIRE(sum3 == sum);
}

TEST_CASE("Separate index spaces", "[newtype]") {
	// A small CSR-style graph: nodes and edge-slots live in different index
	// spaces, and the compiler prevents mixing them up.
	INT_NEWTYPE(node);
	INT_NEWTYPE(edge_slot);

	vec<node, edge_slot> deg(node(3), edge_slot(0));
	std::vector<std::pair<node, node>> edges = {{node(0), node(1)}, {node(1), node(2)}, {node(0), node(2)}};
	for (auto [u, v_] : edges) { deg[u]++; deg[v_]++; }

	vec<node, edge_slot> start(node(4)); // prefix sums of degrees
	for (node u = 0; u < deg.size(); ++u) start[u + 1] = start[u] + deg[u];
	// start[deg.size()] would not compile with a raw int; edge_slot(...) is the edge-space size
	vec<edge_slot, node> adj(start[node(3)]);
	vec<node, edge_slot> cur = start;
	for (auto [u, v_] : edges) {
		adj[cur[u]++] = v_;
		adj[cur[v_]++] = u;
	}
	REQUIRE(adj.size() == edge_slot(6));
	// adj[node(0)] does not compile; only edge_slot keys are accepted
	static_assert(!std::is_invocable_v<decltype([](auto& a2, auto k) -> decltype(a2[k]) { return a2[k]; }), decltype(adj)&, node>);

	int total = 0;
	for (node u = 0; u < deg.size(); ++u) {
		for (edge_slot i = start[u]; i < start[u + 1]; ++i) total += adj[i].v;
	}
	REQUIRE(total == (0 + 1) + (1 + 2) + (0 + 2));
}

TEST_CASE("Values with wider prefix sums", "[newtype]") {
	// Values (not indices) as newtypes: prefix sums may need a wider type.
	INT_NEWTYPE(idx_t);
	INT_NEWTYPE(val_t, int32_t);
	INT_NEWTYPE(sum_t, int64_t);

	vec<idx_t, val_t> a(idx_t(4));
	for (idx_t i = 0; i < a.size(); ++i) a[i] = val_t(1000000000 + i.v);

	vec<idx_t, sum_t> pfx(a.size() + 1);
	for (idx_t i = 0; i < a.size(); ++i) pfx[i + 1] = sum_t(pfx[i].v + a[i].v);
	REQUIRE(pfx.back() == sum_t(int64_t(4000000006)));
}
