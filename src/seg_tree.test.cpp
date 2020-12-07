#include "catch.hpp"

#include "seg_tree.hpp"

#include <type_traits>

TEMPLATE_TEST_CASE("Segment Tree Layouts", "[seg_tree][template]", seg_tree::in_order_layout, seg_tree::circular_layout) {
	for (int N : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 101, 127, 128, 129}) {
		auto seg = TestType(N);
		for (int i = 0; i < N; i++) {
			auto pt = seg.get_point(i);
			REQUIRE(seg.get_leaf_index(pt) == i);
			if constexpr (std::is_same_v<TestType, seg_tree::in_order_layout>) {
				REQUIRE(seg.get_node_bounds(pt) == std::array<int, 2>({i,i+1}));
			} else {
				REQUIRE(seg.get_node_bounds(pt) == std::array<int, 2>({i,(i+1)%N}));
			}
			REQUIRE(seg.get_node_size(pt) == 1);
		}
		for (int a = N-1; a >= 1; a--) {
			auto pt = seg_tree::point(a);
			REQUIRE(seg.get_node_size(pt) == seg.get_node_size(pt.c(0)) + seg.get_node_size(pt.c(1)));
			REQUIRE(seg.get_node_bounds(pt)[0] == seg.get_node_bounds(pt.c(0))[0]);
			REQUIRE(seg.get_node_bounds(pt)[1] == seg.get_node_bounds(pt.c(1))[1]);
			REQUIRE(seg.get_node_bounds(pt.c(0))[1] == seg.get_node_bounds(pt.c(1))[0]);
		}

		for (int l = 0; l <= N; l++) {
			for (int r = l; r <= N; r++) {
				auto rng = seg.get_range(l, r);

				int x = l, y = r;
				if constexpr (std::is_same_v<TestType, seg_tree::circular_layout>) {
					x %= N;
					y %= N;
				}
				rng.for_each([&](auto a, bool d) {
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
		}
	}
}
