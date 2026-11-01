#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <utility>

#include "containers/vec.hpp"
#include "graph/mcmf.hpp"

namespace {

struct edge_spec {
	int from, to, cap;
	int64_t cost;
};

// Successive shortest paths by Bellman-Ford on the residual graph, augmenting
// one unit at a time while the shortest s-t path costs at most max_cost.
std::pair<int, int64_t> reference_mcmf(int N, wala::vec<edge_spec> edges, int s, int t, int64_t max_cost) {
	int E = int(edges.size());
	wala::vec<int> rev_cap(E, 0);
	int tot_flow = 0;
	int64_t tot_cost = 0;
	while (true) {
		const int64_t INF = INT64_MAX / 4;
		wala::vec<int64_t> dist(N, INF);
		wala::vec<int> prv(N, -1);
		dist[s] = 0;
		for (int iter = 0; iter < N; iter++) {
			for (int e = 0; e < E; e++) {
				auto [u, v, cap, cost] = edges[e];
				if (cap > 0 && dist[u] < INF && dist[u] + cost < dist[v]) {
					dist[v] = dist[u] + cost;
					prv[v] = e;
				}
				if (rev_cap[e] > 0 && dist[v] < INF && dist[v] - cost < dist[u]) {
					dist[u] = dist[v] - cost;
					prv[u] = ~e;
				}
			}
		}
		if (dist[t] >= INF || dist[t] > max_cost) break;
		for (int cur = t; cur != s; ) {
			int e = prv[cur];
			if (e >= 0) {
				edges[e].cap--, rev_cap[e]++;
				cur = edges[e].from;
			} else {
				edges[~e].cap++, rev_cap[~e]--;
				cur = edges[~e].to;
			}
		}
		tot_flow++;
		tot_cost += dist[t];
	}
	return {tot_flow, tot_cost};
}

wala::vec<edge_spec> random_graph(std::mt19937& mt, int N, int E) {
	wala::vec<edge_spec> edges(E, edge_spec{});
	for (auto& e : edges) {
		e.from = std::uniform_int_distribution<int>(0, N-1)(mt);
		e.to = std::uniform_int_distribution<int>(0, N-1)(mt);
		e.cap = std::uniform_int_distribution<int>(0, 3)(mt);
		e.cost = std::uniform_int_distribution<int>(0, 4)(mt);
	}
	return edges;
}

template <typename MCMF>
void check_against_reference(std::mt19937& mt) {
	for (int iter = 0; iter < 300; iter++) {
		int N = std::uniform_int_distribution<int>(2, 6)(mt);
		int E = std::uniform_int_distribution<int>(0, 10)(mt);
		auto edges = random_graph(mt, N, E);
		int s = 0, t = N-1;
		int64_t max_cost = std::uniform_int_distribution<int>(0, 2)(mt) ? MCMF::INF_COST - 1 : std::uniform_int_distribution<int>(0, 8)(mt);
		auto expected = reference_mcmf(N, edges, s, t, max_cost);

		MCMF g(N, E);
		for (auto [u, v, cap, cost] : edges) g.add_edge(u, v, cap, cost);
		REQUIRE(g.max_flow(s, t, max_cost) == expected);

		MCMF h(N, E);
		for (auto [u, v, cap, cost] : edges) h.add_edge(u, v, cap, cost);
		auto flows = h.all_flows(s, t, max_cost);
		int tot_flow = 0;
		int64_t tot_cost = 0;
		for (int i = 0; i < int(flows.size()); i++) {
			REQUIRE((i == 0 || flows[i].second >= flows[i-1].second));
			tot_flow += flows[i].first;
			tot_cost += flows[i].first * flows[i].second;
		}
		REQUIRE(std::pair(tot_flow, tot_cost) == expected);
	}
}

// Two calls on one instance (with edges added in between) must agree with a
// single call on an instance that saw all the edges up front: the first
// call's flow is min-cost for its value, and the added edges have nonnegative
// reduced cost, so the second call completes the same min-cost max flow.
template <typename MCMF>
void check_multiple_calls(std::mt19937& mt) {
	for (int iter = 0; iter < 300; iter++) {
		int N = std::uniform_int_distribution<int>(2, 6)(mt);
		int E1 = std::uniform_int_distribution<int>(0, 8)(mt);
		int E2 = std::uniform_int_distribution<int>(0, 4)(mt);
		auto edges = random_graph(mt, N, E1 + E2);
		int s = 0, t = N-1;
		int64_t mid_cost = std::uniform_int_distribution<int>(0, 2)(mt) ? MCMF::INF_COST - 1 : std::uniform_int_distribution<int>(0, 8)(mt);

		MCMF g(N, E1 + E2);
		for (int e = 0; e < E1; e++) g.add_edge(edges[e].from, edges[e].to, edges[e].cap, edges[e].cost);
		auto [f1, c1] = g.max_flow(s, t, mid_cost);
		for (int e = E1; e < E1 + E2; e++) {
			// Potentials from the first call may make the new edge's reduced
			// cost negative; lift its cost to keep add_edge's precondition.
			edges[e].cost = std::max(edges[e].cost, g.pi[edges[e].to] - g.pi[edges[e].from]);
			g.add_edge(edges[e].from, edges[e].to, edges[e].cap, edges[e].cost);
		}
		auto [f2, c2] = g.max_flow(s, t);
		auto [f3, c3] = g.max_flow(s, t);
		REQUIRE(f3 == 0);
		REQUIRE(c3 == 0);

		MCMF h(N, E1 + E2);
		for (auto [u, v, cap, cost] : edges) h.add_edge(u, v, cap, cost);
		REQUIRE(h.max_flow(s, t) == std::pair(f1 + f2, c1 + c2));
	}
}

} // namespace

TEST_CASE("MCMF_SSPA matches reference", "[mcmf]") {
	std::mt19937 mt(48);
	check_against_reference<wala::MCMF_SSPA<int, int64_t>>(mt);
}

TEST_CASE("MCMF_Dinic matches reference", "[mcmf]") {
	std::mt19937 mt(49);
	check_against_reference<wala::MCMF_Dinic<int, int64_t>>(mt);
}

TEST_CASE("MCMF across multiple calls", "[mcmf]") {
	std::mt19937 mt(50);
	check_multiple_calls<wala::MCMF_SSPA<int, int64_t>>(mt);
	check_multiple_calls<wala::MCMF_Dinic<int, int64_t>>(mt);
}

TEST_CASE("Dinic max flow across multiple calls", "[mcmf]") {
	std::mt19937 mt(51);
	for (int iter = 0; iter < 300; iter++) {
		int N = std::uniform_int_distribution<int>(2, 6)(mt);
		int E1 = std::uniform_int_distribution<int>(0, 8)(mt);
		int E2 = std::uniform_int_distribution<int>(0, 4)(mt);
		auto edges = random_graph(mt, N, E1 + E2);
		int s = 0, t = N-1;
		auto expected = reference_mcmf(N, edges, s, t, INT64_MAX / 4).first;

		wala::Dinic<int> g(N, E1 + E2);
		for (int e = 0; e < E1; e++) g.add_edge(edges[e].from, edges[e].to, edges[e].cap);
		int f1 = g.max_flow(s, t);
		for (int e = E1; e < E1 + E2; e++) g.add_edge(edges[e].from, edges[e].to, edges[e].cap);
		int f2 = g.max_flow(s, t);
		REQUIRE(g.max_flow(s, t) == 0);
		REQUIRE(f1 + f2 == expected);
	}
}
