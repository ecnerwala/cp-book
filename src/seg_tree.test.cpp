#include <catch2/catch_template_test_macros.hpp>

#include "seg_tree.hpp"

#include <type_traits>

TEMPLATE_TEST_CASE("Segment Tree Layouts", "[seg_tree][template]", seg_tree::in_order_layout, seg_tree::circular_layout) {
	for (int N : {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 101, 127, 128, 129}) {
		auto seg = TestType(N);
		for (int i = 0; i < N; i++) {
			auto pt = seg.get_point(i);
			REQUIRE(seg.get_leaf_index(pt) == i);
			REQUIRE(seg.get_node_bounds(pt) == std::array<int, 2>({i,i+1}));
			REQUIRE(seg.get_node_size(pt) == 1);
		}
		for (seg_tree::point a(N-1); a >= 1; a--) {
			auto pt = seg_tree::point(a);
			REQUIRE(seg.get_node_size(pt) == seg.get_node_size(pt.c(0)) + seg.get_node_size(pt.c(1)));
			REQUIRE(seg.get_node_bounds(pt)[0] == seg.get_node_bounds(pt.c(0))[0]);
			REQUIRE(seg.get_node_bounds(pt)[1] == seg.get_node_bounds(pt.c(1))[1]);
			if constexpr (std::is_same_v<TestType, seg_tree::in_order_layout>) {
				REQUIRE(seg.get_node_bounds(pt.c(0))[1] == seg.get_node_bounds(pt.c(1))[0]);
			} else {
				REQUIRE(seg.get_node_bounds(pt.c(0))[1] % N == seg.get_node_bounds(pt.c(1))[0]);
			}
		}

		for (int l = 0; l <= N; l++) {
			for (int r = l; r <= N; r++) {
				auto rng = seg.get_range(l, r);

				{
					int x = l, y = r;
					rng.for_each([&](auto a) {
						auto bounds = seg.get_node_bounds(a);
						if (x == bounds[0]) {
							x = bounds[1];
						} else if (y == bounds[1]) {
							y = bounds[0];
						} else assert(false);
					});
					REQUIRE(x == y);
				}
				{
					int x = l, y = r;
					rng.for_each_with_side([&](auto a, bool d) {
						auto bounds = seg.get_node_bounds(a);
						if (d == 0) {
							REQUIRE(x == bounds[0]);
							x = bounds[1];
						} else if (d == 1) {
							REQUIRE(y == bounds[1]);
							y = bounds[0];
						} else assert(false);
					});
					REQUIRE(x == y);
				}
				{
					int x = l;
					rng.for_each_l_to_r([&](auto a) {
						auto bounds = seg.get_node_bounds(a);
						REQUIRE(x == bounds[0]);
						x = bounds[1];
					});
					REQUIRE(x == r);
				}
				{
					int y = r;
					rng.for_each_r_to_l([&](auto a) {
						auto bounds = seg.get_node_bounds(a);
						REQUIRE(y == bounds[1]);
						y = bounds[0];
					});
					REQUIRE(y == l);
				}
			}
		}
	}
}
