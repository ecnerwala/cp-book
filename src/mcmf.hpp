#pragma once
#include<bits/stdc++.h>
#include<bits/extc++.h>

// NOTE: This doesn't support negative-cost edges; you can adjust edge weights
// (e.g. by precomputing a potential function) to make them positive.

template <typename flow_t = int, typename cost_t = int64_t>
struct MCMF_SSPA {
	int N;
	std::vector<std::vector<int>> adj;
	struct edge_t {
		int dest;
		flow_t cap;
		cost_t cost;
	};
	std::vector<edge_t> edges;

	std::vector<char> seen;
	std::vector<cost_t> pi;
	std::vector<int> prv;

	explicit MCMF_SSPA(int N_) : N(N_), adj(N), pi(N, 0), prv(N) {}

	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		int e = int(edges.size());
		edges.emplace_back(edge_t{to, cap, cost});
		edges.emplace_back(edge_t{from, 0, -cost});
		adj[from].push_back(e);
		adj[to].push_back(e+1);
	}

	const cost_t INF_COST = std::numeric_limits<cost_t>::max() / 4;
	const flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 4;
	std::vector<cost_t> dist;
	__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
	std::vector<typename decltype(q)::point_iterator> its;
	cost_t dijkstra(int s, int t) {
		dist.assign(N, INF_COST);
		dist[s] = 0;

		its.assign(N, q.end());
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

		swap(pi, dist);
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

	std::pair<flow_t, cost_t> max_flow(int s, int t) {
		assert(s != t);
		flow_t tot_flow = 0; cost_t tot_cost = 0;
		while (dijkstra(s, t) < INF_COST) {
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
	std::vector<std::vector<int>> adj;
	struct edge_t {
		int dest;
		flow_t cap;
		cost_t cost;
	};
	std::vector<edge_t> edges;

	std::vector<char> seen;
	std::vector<cost_t> pi;

	explicit MCMF_Dinic(int N_) : N(N_), adj(N), pi(N, 0) {}

	void add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		int e = int(edges.size());
		edges.emplace_back(edge_t{to, cap, cost});
		edges.emplace_back(edge_t{from, 0, -cost});
		adj[from].push_back(e);
		adj[to].push_back(e+1);
	}

	static constexpr cost_t INF_COST = std::numeric_limits<cost_t>::max() / 4;
	static constexpr flow_t INF_FLOW = std::numeric_limits<flow_t>::max() / 4;
	std::vector<cost_t> dist;
	__gnu_pbds::priority_queue<std::pair<cost_t, int>> q;
	std::vector<typename decltype(q)::point_iterator> its;
	cost_t dijkstra(int s, int t) {
		dist.assign(N, INF_COST);
		dist[s] = 0;

		its.assign(N, q.end());
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

	std::vector<int> buf;
	std::vector<int> level;
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
			buf.reserve(N);
			level.assign(N, -1);
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
			buf.assign(N, 0);
			tot_flow += dinic_dfs(s, t, INF_FLOW);
		}
		return tot_flow;
	}

	std::vector<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		std::vector<std::pair<flow_t, cost_t>> res;
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
			int cur_flow = dinic(s, t);
			tot_flow += cur_flow;
			tot_cost += cur_flow * pi[t];
		}
		return {tot_flow, tot_cost};
	}
};
