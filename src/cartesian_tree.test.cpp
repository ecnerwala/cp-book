#include "cartesian_tree.hpp"

#include "catch.hpp"
#include <bits/stdc++.h>

TEST_CASE("Cartesian Tree", "[cartesian_tree]") {
	std::mt19937 mt(48);
	for (int sz : {0, 1, 2, 3, 5, 8, 13}) {
		std::vector<int> v(sz);
		iota(v.begin(), v.end(), 0);
		shuffle(v.begin(), v.end(), mt);
		{
			CartesianTree t = CartesianTree::build_min_tree(v);
			for (int i = 1; i < int(t.nodes.size()); i += 2) {
				auto cur = &t.nodes[i];
				REQUIRE(cur->m == i/2);
				REQUIRE(cur->l <= cur->m);
				REQUIRE(cur->m <= cur->r);

				REQUIRE(cur->c[0]->l == cur->l);
				REQUIRE(cur->c[0]->r == cur->m-1);

				REQUIRE(cur->c[1]->l == cur->m+1);
				REQUIRE(cur->c[1]->r == cur->r);

				REQUIRE((cur->c[0]->l > cur->c[0]->r || v[cur->m] < v[cur->c[0]->m]));
				REQUIRE((cur->c[1]->l > cur->c[1]->r || v[cur->m] < v[cur->c[1]->m]));
			}
		}
		{
			CartesianTree t = CartesianTree::build_max_tree(v);
			for (int i = 1; i < int(t.nodes.size()); i += 2) {
				auto cur = &t.nodes[i];
				REQUIRE(cur->m == i/2);
				REQUIRE(cur->l <= cur->m);
				REQUIRE(cur->m <= cur->r);

				REQUIRE(cur->c[0]->l == cur->l);
				REQUIRE(cur->c[0]->r == cur->m-1);

				REQUIRE(cur->c[1]->l == cur->m+1);
				REQUIRE(cur->c[1]->r == cur->r);

				REQUIRE((cur->c[0]->l > cur->c[0]->r || v[cur->m] > v[cur->c[0]->m]));
				REQUIRE((cur->c[1]->l > cur->c[1]->r || v[cur->m] > v[cur->c[1]->m]));
			}
		}
	}
}
