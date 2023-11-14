#include "cartesian_tree.hpp"

#include <catch2/catch_test_macros.hpp>
#include <bits/stdc++.h>

TEST_CASE("Cartesian Tree", "[cartesian_tree]") {
	std::mt19937 mt(48);
	for (int sz : {0, 1, 2, 3, 5, 8, 13}) {
		std::vector<int> v(sz);
		iota(v.begin(), v.end(), 0);
		shuffle(v.begin(), v.end(), mt);
		{
			CartesianTree t = CartesianTree::build_min_tree(v);
			for (int i = 1; i < int(t.size()); i += 2) {
				REQUIRE(t[i].m == i/2);
				REQUIRE(t[i].l <= t[i].m);
				REQUIRE(t[i].m <= t[i].r);

				REQUIRE(t[t[i].c[0]].l == t[i].l);
				REQUIRE(t[t[i].c[0]].r == t[i].m-1);

				REQUIRE(t[t[i].c[1]].l == t[i].m+1);
				REQUIRE(t[t[i].c[1]].r == t[i].r);

				REQUIRE((t[t[i].c[0]].l > t[t[i].c[0]].r || v[t[i].m] < v[t[t[i].c[0]].m]));
				REQUIRE((t[t[i].c[1]].l > t[t[i].c[1]].r || v[t[i].m] < v[t[t[i].c[1]].m]));
			}
		}
		{
			CartesianTree t = CartesianTree::build_max_tree(v);
			for (int i = 1; i < int(t.size()); i += 2) {
				REQUIRE(t[i].m == i/2);
				REQUIRE(t[i].l <= t[i].m);
				REQUIRE(t[i].m <= t[i].r);

				REQUIRE(t[t[i].c[0]].l == t[i].l);
				REQUIRE(t[t[i].c[0]].r == t[i].m-1);

				REQUIRE(t[t[i].c[1]].l == t[i].m+1);
				REQUIRE(t[t[i].c[1]].r == t[i].r);

				REQUIRE((t[t[i].c[0]].l > t[t[i].c[0]].r || v[t[i].m] > v[t[t[i].c[0]].m]));
				REQUIRE((t[t[i].c[1]].l > t[t[i].c[1]].r || v[t[i].m] > v[t[t[i].c[1]].m]));
			}
		}
	}
}
