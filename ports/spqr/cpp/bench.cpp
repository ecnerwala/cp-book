#include "graph/spqr_tree.hpp"
#include <chrono>
#include <cstdio>
#include <iostream>
int main(int argc, char** argv) {
	int R = argc > 1 ? atoi(argv[1]) : 5;
	int NV, NE; std::cin >> NV >> NE;
	std::vector<std::array<int, 2>> edges(NE);
	for (auto& [u, v] : edges) std::cin >> u >> v;
	auto bench = [&](const char* name, auto f) {
		double best = 1e18; long long sink = 0;
		for (int r = 0; r < R; r++) {
			auto t0 = std::chrono::steady_clock::now();
			auto t = f();
			auto t1 = std::chrono::steady_clock::now();
			sink += t.par.size();
			best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
		}
		printf("cpp  %-10s %8.2f ms  (items=%lld)\n", name, best, sink / R);
	};
	bench("spqr", [&] { return wala::spqr_tree::build(NV, edges); });
	bench("planar", [&] { return wala::planar_spqr_tree::build(NV, edges); });
}
