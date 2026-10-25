#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <utility>
#include <ext/pb_ds/priority_queue.hpp>

#include "containers/csr.hpp"
#include "containers/vec.hpp"

namespace wala {

// NOTE: This doesn't support negative-cost edges; you can adjust edge weights
// (e.g. by precomputing a potential function) to make them positive.
//
// Usage: construct with (N, max_edges) and add_edge up to max_edges times
// (total, across all calls).
// Each max_flow / all_flows builds its adjacency from the current edge list
// and continues from the current residual capacities and potentials, so
// calls can be interleaved with add_edge.

namespace detail {
template <typename edge_t>
csr<int> flow_adjacency(int N, const bounded_vec<edge_t>& edges) {
	csr_index_builder idx_builder(N);
	for (int e = 0; e < int(edges.size()); e++) idx_builder.count(edges[e^1].dest);
	csr_builder<int> builder(std::move(idx_builder));
	for (int e = 0; e < int(edges.size()); e++) builder.push(edges[e^1].dest) = e;
	return std::move(builder).finalize();
}
} // namespace detail

template <typename flow_t = int, typename cost_t = int64_t>
struct MCMF_SSPA {
	int N;
	struct edge_t {
		int dest;
		flow_t cap;
		cost_t cost;
	};
	bounded_vec<edge_t> edges;

	vec<cost_t> pi;

	explicit MCMF_SSPA(int N_, int max_edges) : N(N_), edges(with_capacity{2 * max_edges}), pi(N, 0) {}

	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		edges.push_back(edge_t{to, cap, cost});
		edges.push_back(edge_t{from, 0, -cost});
	}

	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() / 4;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 4;

	struct solver {
		MCMF_SSPA& g;
		csr<int> adj;
		vec<cost_t> dist;
		vec<int> prv;
		__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
		vec<typename decltype(q)::point_iterator> its;

		explicit solver(MCMF_SSPA& g_) : g(g_), adj(detail::flow_adjacency(g.N, g.edges)), dist(g.N), prv(g.N), its(g.N) {}

		cost_t dijkstra(int s, int t) {
			std::ranges::fill(dist, INF_COST);
			dist[s] = 0;

			std::ranges::fill(its, q.end());
			its[s] = q.push({-(dist[s] - g.pi[s]), s});

			while (!q.empty()) {
				int i = q.top().second; q.pop();
				cost_t d = dist[i];
				for (int e : adj[i]) {
					if (g.edges[e].cap) {
						int j = g.edges[e].dest;
						cost_t nd = d + g.edges[e].cost;
						if (nd < dist[j]) {
							dist[j] = nd;
							prv[j] = e;
							if (its[j] == q.end()) {
								its[j] = q.push({-(dist[j] - g.pi[j]), j});
							} else {
								q.modify(its[j], {-(dist[j] - g.pi[j]), j});
							}
						}
					}
				}
			}

			std::swap(g.pi, dist);
			return g.pi[t];
		}

		flow_t path(int s, int t) {
			flow_t cur_flow = std::numeric_limits<flow_t>::max();
			for (int cur = t; cur != s; ) {
				int e = prv[cur];
				int nxt = g.edges[e^1].dest;
				cur_flow = std::min(cur_flow, g.edges[e].cap);
				cur = nxt;
			}
			for (int cur = t; cur != s; ) {
				int e = prv[cur];
				int nxt = g.edges[e^1].dest;
				g.edges[e].cap -= cur_flow;
				g.edges[e^1].cap += cur_flow;
				cur = nxt;
			}
			return cur_flow;
		}
	};

	growable_vec<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		solver run(*this);
		growable_vec<std::pair<flow_t, cost_t>> res;
		while (run.dijkstra(s, t) <= max_cost) {
			assert(res.empty() || pi[t] >= res.back().second);
			flow_t f = run.path(s, t);
			res.push_back({f, pi[t]});
		}
		return res;
	}

	std::pair<flow_t, cost_t> max_flow(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		solver run(*this);
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		while (run.dijkstra(s, t) <= max_cost) {
			flow_t cur_flow = run.path(s, t);
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
		}
		return {tot_flow, tot_cost};
	}
};

template <typename flow_t = int, typename cost_t = int64_t>
struct MCMF_Dinic {
	int N;
	struct edge_t {
		int dest;
		flow_t cap;
		cost_t cost;
	};
	bounded_vec<edge_t> edges;

	vec<cost_t> pi;

	explicit MCMF_Dinic(int N_, int max_edges) : N(N_), edges(with_capacity{2 * max_edges}), pi(N, 0) {}

	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		edges.push_back(edge_t{to, cap, cost});
		edges.push_back(edge_t{from, 0, -cost});
	}

	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() / 4;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 4;

	struct solver {
		MCMF_Dinic& g;
		csr<int> adj;
		vec<cost_t> dist;
		__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
		vec<typename decltype(q)::point_iterator> its;
		bounded_vec<int> buf;
		vec<int> level;

		explicit solver(MCMF_Dinic& g_) : g(g_), adj(detail::flow_adjacency(g.N, g.edges)), dist(g.N), its(g.N), buf(with_capacity{g.N}), level(g.N) {}

		cost_t dijkstra(int s, int t) {
			std::ranges::fill(dist, INF_COST);
			dist[s] = 0;

			std::ranges::fill(its, q.end());
			its[s] = q.push({-(dist[s] - g.pi[s]), s});

			while (!q.empty()) {
				int i = q.top().second; q.pop();
				cost_t d = dist[i];
				for (int e : adj[i]) {
					if (g.edges[e].cap) {
						int j = g.edges[e].dest;
						cost_t nd = d + g.edges[e].cost;
						if (nd < dist[j]) {
							dist[j] = nd;
							if (its[j] == q.end()) {
								its[j] = q.push({-(dist[j] - g.pi[j]), j});
							} else {
								q.modify(its[j], {-(dist[j] - g.pi[j]), j});
							}
						}
					}
				}
			}

			std::swap(g.pi, dist);
			return g.pi[t];
		}

		flow_t dinic_dfs(int cur, int t, flow_t f) {
			if (cur == t) return f;
			flow_t cur_f = 0;
			assert(f > 0);
			for (; buf[cur] < int(adj[cur].size()); buf[cur]++) {
				int e = adj[cur][buf[cur]];
				int nxt = g.edges[e].dest;
				if (level[nxt] == level[cur] + 1 && g.edges[e].cap > 0 && g.edges[e].cost == g.pi[nxt] - g.pi[cur]) {
					flow_t v = dinic_dfs(nxt, t, std::min(f, g.edges[e].cap));
					g.edges[e].cap -= v;
					g.edges[e^1].cap += v;
					f -= v;
					cur_f += v;
					if (f == 0) break;
				}
			}
			return cur_f;
		}
		flow_t dinic(int s, int t) {
			flow_t tot_flow = 0;
			while (true) {
				buf.clear();
				std::ranges::fill(level, -1);
				buf.push_back(s);
				level[s] = 0;
				for (int z = 0; z < int(buf.size()); z++) {
					int cur = buf[z];
					for (int e : adj[cur]) {
						int nxt = g.edges[e].dest;
						if (g.edges[e].cap > 0 && g.edges[e].cost == g.pi[nxt] - g.pi[cur] && level[nxt] == -1) {
							level[nxt] = level[cur] + 1;
							buf.push_back(nxt);
						}
					}
				}
				if (level[t] == -1) break;
				buf.clear();
				buf.grow_to(g.N, 0);
				tot_flow += dinic_dfs(s, t, INF_FLOW);
			}
			return tot_flow;
		}
	};

	growable_vec<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		solver run(*this);
		growable_vec<std::pair<flow_t, cost_t>> res;
		while (run.dijkstra(s, t) <= max_cost) {
			assert(res.empty() || pi[t] > res.back().second);
			flow_t f = run.dinic(s, t);
			res.push_back({f, pi[t]});
		}
		return res;
	}

	std::pair<flow_t, cost_t> max_flow(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		solver run(*this);
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		while (run.dijkstra(s, t) <= max_cost) {
			flow_t cur_flow = run.dinic(s, t);
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
		}
		return {tot_flow, tot_cost};
	}
};

template <typename flow_t = int, typename tot_flow_t = flow_t>
struct Dinic {
	int N;
	struct edge_t {
		int dest;
		flow_t cap;
	};
	bounded_vec<edge_t> edges;

	explicit Dinic(int N_, int max_edges) : N(N_), edges(with_capacity{2 * max_edges}) {}

	void add_edge(int from, int to, flow_t cap) {
		return add_bi_edge(from, to, cap, 0);
	}

	void add_bi_edge(int from, int to, flow_t cap, flow_t rev_cap) {
		assert(cap >= 0);
		assert(rev_cap >= 0);
		edges.push_back(edge_t{to, cap});
		edges.push_back(edge_t{from, rev_cap});
	}

	static constexpr tot_flow_t INF_FLOW = std::numeric_limits<tot_flow_t>::max() / 4;

	struct solver {
		Dinic& g;
		csr<int> adj;
		bounded_vec<int> buf;
		vec<int> level;

		explicit solver(Dinic& g_) : g(g_), adj(detail::flow_adjacency(g.N, g.edges)), buf(with_capacity{g.N}), level(g.N) {}

		tot_flow_t dinic_dfs(int cur, int t, tot_flow_t f) {
			if (cur == t) return f;
			tot_flow_t cur_f = 0;
			assert(f > 0);
			for (; buf[cur] < int(adj[cur].size()); buf[cur]++) {
				int e = adj[cur][buf[cur]];
				int nxt = g.edges[e].dest;
				if (level[nxt] == level[cur] + 1 && g.edges[e].cap > 0) {
					flow_t v = flow_t(dinic_dfs(nxt, t, std::min<tot_flow_t>(f, g.edges[e].cap)));
					g.edges[e].cap -= v;
					g.edges[e^1].cap += v;
					f -= v;
					cur_f += v;
					if (f == 0) break;
				}
			}
			return cur_f;
		}
		tot_flow_t dinic(int s, int t) {
			tot_flow_t tot_flow = 0;
			while (true) {
				buf.clear();
				std::ranges::fill(level, -1);
				buf.push_back(s);
				level[s] = 0;
				for (int z = 0; z < int(buf.size()); z++) {
					int cur = buf[z];
					for (int e : adj[cur]) {
						int nxt = g.edges[e].dest;
						if (g.edges[e].cap > 0 && level[nxt] == -1) {
							level[nxt] = level[cur] + 1;
							buf.push_back(nxt);
						}
					}
				}
				if (level[t] == -1) break;
				buf.clear();
				buf.grow_to(g.N, 0);
				tot_flow += dinic_dfs(s, t, INF_FLOW);
			}
			return tot_flow;
		}
	};

	tot_flow_t max_flow(int s, int t) {
		solver run(*this);
		return run.dinic(s, t);
	}
};

} // namespace wala
