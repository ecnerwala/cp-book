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
// Usage: construct with (N, max_edges), add_edge up to max_edges times, then
// build() once before any query.

template <typename flow_t = int, typename cost_t = int64_t>
struct MCMF_SSPA {
	int N;
	struct edge_t {
		int dest;
		flow_t cap;
		cost_t cost;
	};
	bounded_vec<edge_t> edges;
	csr_index_builder adj_builder;
	csr<int> adj;

	vec<cost_t> pi;
	vec<int> prv;

	explicit MCMF_SSPA(int N_, int max_edges) : N(N_), edges(with_capacity{2 * max_edges}), adj_builder(N), pi(N, 0), prv(N), dist(N), its(N) {}

	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		edges.push_back(edge_t{to, cap, cost});
		edges.push_back(edge_t{from, 0, -cost});
		adj_builder.count(from);
		adj_builder.count(to);
	}

	void build() {
		csr_builder<int> builder(std::move(adj_builder));
		for (int e = 0; e < int(edges.size()); e++) builder.push(edges[e^1].dest) = e;
		adj = std::move(builder).finalize();
	}

	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() / 4;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 4;
	vec<cost_t> dist;
	__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
	vec<typename decltype(q)::point_iterator> its;
	cost_t dijkstra(int s, int t) {
		std::ranges::fill(dist, INF_COST);
		dist[s] = 0;

		std::ranges::fill(its, q.end());
		its[s] = q.push({-(dist[s] - pi[s]), s});

		while (!q.empty()) {
			int i = q.top().second; q.pop();
			cost_t d = dist[i];
			for (int e : adj[i]) {
				if (edges[e].cap) {
					int j = edges[e].dest;
					cost_t nd = d + edges[e].cost;
					if (nd < dist[j]) {
						dist[j] = nd;
						prv[j] = e;
						if (its[j] == q.end()) {
							its[j] = q.push({-(dist[j] - pi[j]), j});
						} else {
							q.modify(its[j], {-(dist[j] - pi[j]), j});
						}
					}
				}
			}
		}

		std::swap(pi, dist);
		return pi[t];
	}

	flow_t path(int s, int t) {
		flow_t cur_flow = std::numeric_limits<flow_t>::max();
		for (int cur = t; cur != s; ) {
			int e = prv[cur];
			int nxt = edges[e^1].dest;
			cur_flow = std::min(cur_flow, edges[e].cap);
			cur = nxt;
		}
		for (int cur = t; cur != s; ) {
			int e = prv[cur];
			int nxt = edges[e^1].dest;
			edges[e].cap -= cur_flow;
			edges[e^1].cap += cur_flow;
			cur = nxt;
		}
		return cur_flow;
	}

	growable_vec<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		growable_vec<std::pair<flow_t, cost_t>> res;
		while (dijkstra(s, t) <= max_cost) {
			assert(res.empty() || pi[t] >= res.back().second);
			flow_t f = path(s, t);
			res.push_back({f, pi[t]});
		}
		return res;
	}

	std::pair<flow_t, cost_t> max_flow(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		while (dijkstra(s, t) <= max_cost) {
			flow_t cur_flow = path(s, t);
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
	csr_index_builder adj_builder;
	csr<int> adj;

	vec<cost_t> pi;

	explicit MCMF_Dinic(int N_, int max_edges) : N(N_), edges(with_capacity{2 * max_edges}), adj_builder(N), pi(N, 0), dist(N), its(N), buf(with_capacity{N}), level(N) {}

	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		edges.push_back(edge_t{to, cap, cost});
		edges.push_back(edge_t{from, 0, -cost});
		adj_builder.count(from);
		adj_builder.count(to);
	}

	void build() {
		csr_builder<int> builder(std::move(adj_builder));
		for (int e = 0; e < int(edges.size()); e++) builder.push(edges[e^1].dest) = e;
		adj = std::move(builder).finalize();
	}

	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() / 4;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 4;
	vec<cost_t> dist;
	__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
	vec<typename decltype(q)::point_iterator> its;
	cost_t dijkstra(int s, int t) {
		std::ranges::fill(dist, INF_COST);
		dist[s] = 0;

		std::ranges::fill(its, q.end());
		its[s] = q.push({-(dist[s] - pi[s]), s});

		while (!q.empty()) {
			int i = q.top().second; q.pop();
			cost_t d = dist[i];
			for (int e : adj[i]) {
				if (edges[e].cap) {
					int j = edges[e].dest;
					cost_t nd = d + edges[e].cost;
					if (nd < dist[j]) {
						dist[j] = nd;
						if (its[j] == q.end()) {
							its[j] = q.push({-(dist[j] - pi[j]), j});
						} else {
							q.modify(its[j], {-(dist[j] - pi[j]), j});
						}
					}
				}
			}
		}

		std::swap(pi, dist);
		return pi[t];
	}

	bounded_vec<int> buf;
	vec<int> level;
	flow_t dinic_dfs(int cur, int t, flow_t f) {
		if (cur == t) return f;
		flow_t cur_f = 0;
		assert(f > 0);
		for (; buf[cur] < int(adj[cur].size()); buf[cur]++) {
			int e = adj[cur][buf[cur]];
			int nxt = edges[e].dest;
			if (level[nxt] == level[cur] + 1 && edges[e].cap > 0 && edges[e].cost == pi[nxt] - pi[cur]) {
				flow_t v = dinic_dfs(nxt, t, std::min(f, edges[e].cap));
				edges[e].cap -= v;
				edges[e^1].cap += v;
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
					int nxt = edges[e].dest;
					if (edges[e].cap > 0 && edges[e].cost == pi[nxt] - pi[cur] && level[nxt] == -1) {
						level[nxt] = level[cur] + 1;
						buf.push_back(nxt);
					}
				}
			}
			if (level[t] == -1) break;
			buf.clear();
			buf.grow_to(N, 0);
			tot_flow += dinic_dfs(s, t, INF_FLOW);
		}
		return tot_flow;
	}

	growable_vec<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		growable_vec<std::pair<flow_t, cost_t>> res;
		while (dijkstra(s, t) <= max_cost) {
			assert(res.empty() || pi[t] > res.back().second);
			flow_t f = dinic(s, t);
			res.push_back({f, pi[t]});
		}
		return res;
	}

	std::pair<flow_t, cost_t> max_flow(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		while (dijkstra(s, t) <= max_cost) {
			flow_t cur_flow = dinic(s, t);
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
	csr_index_builder adj_builder;
	csr<int> adj;

	explicit Dinic(int N_, int max_edges) : N(N_), edges(with_capacity{2 * max_edges}), adj_builder(N), buf(with_capacity{N}), level(N) {}

	void add_edge(int from, int to, flow_t cap) {
		return add_bi_edge(from, to, cap, 0);
	}

	void add_bi_edge(int from, int to, flow_t cap, flow_t rev_cap) {
		assert(cap >= 0);
		assert(rev_cap >= 0);
		edges.push_back(edge_t{to, cap});
		edges.push_back(edge_t{from, rev_cap});
		adj_builder.count(from);
		adj_builder.count(to);
	}

	void build() {
		csr_builder<int> builder(std::move(adj_builder));
		for (int e = 0; e < int(edges.size()); e++) builder.push(edges[e^1].dest) = e;
		adj = std::move(builder).finalize();
	}

	static constexpr tot_flow_t INF_FLOW = std::numeric_limits<tot_flow_t>::max() / 4;
	bounded_vec<int> buf;
	vec<int> level;
	tot_flow_t dinic_dfs(int cur, int t, tot_flow_t f) {
		if (cur == t) return f;
		tot_flow_t cur_f = 0;
		assert(f > 0);
		for (; buf[cur] < int(adj[cur].size()); buf[cur]++) {
			int e = adj[cur][buf[cur]];
			int nxt = edges[e].dest;
			if (level[nxt] == level[cur] + 1 && edges[e].cap > 0) {
				flow_t v = flow_t(dinic_dfs(nxt, t, std::min<tot_flow_t>(f, edges[e].cap)));
				edges[e].cap -= v;
				edges[e^1].cap += v;
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
					int nxt = edges[e].dest;
					if (edges[e].cap > 0 && level[nxt] == -1) {
						level[nxt] = level[cur] + 1;
						buf.push_back(nxt);
					}
				}
			}
			if (level[t] == -1) break;
			buf.clear();
			buf.grow_to(N, 0);
			tot_flow += dinic_dfs(s, t, INF_FLOW);
		}
		return tot_flow;
	}
	tot_flow_t max_flow(int s, int t) { return dinic(s, t); }
};

} // namespace wala
