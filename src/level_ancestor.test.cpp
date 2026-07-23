#include <catch2/catch_test_macros.hpp>

#include "level_ancestor.hpp"

#include <vector>

namespace {

INT_NEWTYPE(node_t);

template <typename LA, typename Node> void check_tree(const LA& la) {
	using pre_t = typename LA::pre_t;
	//        0
	//       / \
	//      1   2
	//     / \
	//    3   4
	//    |
	//    5
	REQUIRE(la.get_ancestor(Node(5), 0) == Node(5));
	REQUIRE(la.get_ancestor(Node(5), 1) == Node(3));
	REQUIRE(la.get_ancestor(Node(5), 2) == Node(1));
	REQUIRE(la.get_ancestor(Node(5), 3) == Node(0));
	REQUIRE(la.get_ancestor(Node(5), 4) == Node(-1));
	REQUIRE(la.get_ancestor(Node(4), 1) == Node(1));
	REQUIRE(la.get_ancestor(Node(2), 1) == Node(0));

	REQUIRE(la.lca(Node(3), Node(4)) == Node(1));
	REQUIRE(la.lca(Node(5), Node(4)) == Node(1));
	REQUIRE(la.lca(Node(4), Node(2)) == Node(0));
	REQUIRE(la.lca(Node(5), Node(5)) == Node(5));
	REQUIRE(la.lca(Node(0), Node(5)) == Node(0));

	REQUIRE(la.dist(Node(3), Node(4)) == pre_t(2));
	REQUIRE(la.dist(Node(5), Node(2)) == pre_t(4));
	REQUIRE(la.dist(Node(5), Node(1)) == pre_t(2));
	REQUIRE(la.dist(Node(0), Node(0)) == pre_t(0));
}

} // namespace

TEST_CASE("Level ancestor with plain int nodes", "[level_ancestor]") {
	std::vector<int> par = {-1, 0, 0, 1, 1, 3};
	ecnerwala::level_ancestor<> la(par);
	check_tree<decltype(la), int>(la);
}

TEST_CASE("Level ancestor with newtype nodes", "[level_ancestor]") {
	vec<node_t, node_t> par{node_t(-1), node_t(0), node_t(0), node_t(1), node_t(1), node_t(3)};
	ecnerwala::level_ancestor<node_t> la(par);
	check_tree<decltype(la), node_t>(la);
	// la.get_ancestor(0, 1) with a raw runtime int does not compile; literals do:
	REQUIRE(la.get_ancestor(2, 1) == node_t(0));
}
