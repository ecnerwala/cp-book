#pragma once
#include<bits/stdc++.h>
// #include<bits/extc++.h>
#include <ext/pb_ds/priority_queue.hpp>

// NOTE: This doesn't support negative-cost edges; you can adjust edge weights
// (e.g. by precomputing a potential function) to make them positive.

// Edges are staged by add_edge and packed into a CSR adjacency on the first
// solve: each node's edge records (dest, rev, cap, ...) are stored inline and
// contiguously, so the search loops scan flat memory with no indirection
// through a separate adjacency list.
// add_edge returns an edge id; flow(id) reads back the flow after solving.
// No edges may be added once a solve has started.

template <typename flow_t = int, typename cost_t = int64_t>
struct MCMF_SSPA {
	int N;
	struct edge_t {
		int dest;
		int rev;
		flow_t cap;
		cost_t cost;
	};
	struct input_edge_t {
		int from;
		int to;
		flow_t cap;
		cost_t cost;
	};
	std::vector<input_edge_t> input;
	std::vector<int> adj_start;
	std::vector<edge_t> csr;
	std::vector<int> edge_pos; // csr index of each input edge

	std::vector<cost_t> pi;
	std::vector<int> prv;

	explicit MCMF_SSPA(int N_) : N(N_), pi(N, 0), prv(N) {}

	int add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(adj_start.empty());
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		input.push_back(input_edge_t{from, to, cap, cost});
		return int(input.size()) - 1;
	}

	void build() {
		if (!adj_start.empty()) return;
		adj_start.assign(N+1, 0);
		for (const auto& e : input) {
			adj_start[e.from+1]++;
			adj_start[e.to+1]++;
		}
		for (int i = 0; i < N; i++) adj_start[i+1] += adj_start[i];
		csr.resize(2 * input.size());
		edge_pos.resize(input.size());
		std::vector<int> pos(adj_start.begin(), adj_start.end()-1);
		for (int i = 0; i < int(input.size()); i++) {
			const auto& e = input[i];
			int ef = pos[e.from]++, et = pos[e.to]++;
			csr[ef] = edge_t{e.to, et, e.cap, e.cost};
			csr[et] = edge_t{e.from, ef, 0, -e.cost};
			edge_pos[i] = ef;
		}
	}

	flow_t flow(int e) const {
		return input[e].cap - csr[edge_pos[e]].cap;
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
			for (int a = adj_start[i]; a < adj_start[i+1]; a++) {
				const edge_t& e = csr[a];
				if (e.cap) {
					int j = e.dest;
					cost_t nd = d + e.cost;
					if (nd < dist[j]) {
						dist[j] = nd;
						prv[j] = a;
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
			int nxt = csr[csr[e].rev].dest;
			cur_flow = std::min(cur_flow, csr[e].cap);
			cur = nxt;
		}
		for (int cur = t; cur != s; ) {
			int e = prv[cur];
			int nxt = csr[csr[e].rev].dest;
			csr[e].cap -= cur_flow;
			csr[csr[e].rev].cap += cur_flow;
			cur = nxt;
		}
		return cur_flow;
	}

	std::vector<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		build();
		std::vector<std::pair<flow_t, cost_t>> res;
		while (dijkstra(s, t) <= max_cost) {
			assert(res.empty() || pi[t] >= res.back().second);
			flow_t f = path(s, t);
			res.push_back({f, pi[t]});
		}
		return res;
	}

	std::pair<flow_t, cost_t> max_flow(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		build();
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
		int rev;
		flow_t cap;
		cost_t cost;
	};
	struct input_edge_t {
		int from;
		int to;
		flow_t cap;
		cost_t cost;
	};
	std::vector<input_edge_t> input;
	std::vector<int> adj_start;
	std::vector<edge_t> csr;
	std::vector<int> edge_pos;

	std::vector<cost_t> pi;

	explicit MCMF_Dinic(int N_) : N(N_), pi(N, 0) {}

	int add_edge(int from, int to, flow_t cap, cost_t cost) {
		assert(adj_start.empty());
		assert(cap >= 0);
		assert(cost + pi[from] - pi[to] >= 0); // TODO: Remove this restriction
		input.push_back(input_edge_t{from, to, cap, cost});
		return int(input.size()) - 1;
	}

	void build() {
		if (!adj_start.empty()) return;
		adj_start.assign(N+1, 0);
		for (const auto& e : input) {
			adj_start[e.from+1]++;
			adj_start[e.to+1]++;
		}
		for (int i = 0; i < N; i++) adj_start[i+1] += adj_start[i];
		csr.resize(2 * input.size());
		edge_pos.resize(input.size());
		std::vector<int> pos(adj_start.begin(), adj_start.end()-1);
		for (int i = 0; i < int(input.size()); i++) {
			const auto& e = input[i];
			int ef = pos[e.from]++, et = pos[e.to]++;
			csr[ef] = edge_t{e.to, et, e.cap, e.cost};
			csr[et] = edge_t{e.from, ef, 0, -e.cost};
			edge_pos[i] = ef;
		}
	}

	flow_t flow(int e) const {
		return input[e].cap - csr[edge_pos[e]].cap;
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
			for (int a = adj_start[i]; a < adj_start[i+1]; a++) {
				const edge_t& e = csr[a];
				if (e.cap) {
					int j = e.dest;
					cost_t nd = d + e.cost;
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
	std::vector<int> it;
	flow_t dinic_dfs(int cur, int t, flow_t f) {
		if (cur == t) return f;
		flow_t cur_f = 0;
		assert(f > 0);
		for (; it[cur] < adj_start[cur+1]; it[cur]++) {
			edge_t& e = csr[it[cur]];
			int nxt = e.dest;
			if (level[nxt] == level[cur] + 1 && e.cap > 0 && e.cost == pi[nxt] - pi[cur]) {
				flow_t v = dinic_dfs(nxt, t, std::min(f, e.cap));
				e.cap -= v;
				csr[e.rev].cap += v;
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
				for (int a = adj_start[cur]; a < adj_start[cur+1]; a++) {
					const edge_t& e = csr[a];
					int nxt = e.dest;
					if (e.cap > 0 && e.cost == pi[nxt] - pi[cur] && level[nxt] == -1) {
						level[nxt] = level[cur] + 1;
						buf.push_back(nxt);
					}
				}
			}
			if (level[t] == -1) break;
			it.assign(adj_start.begin(), adj_start.end()-1);
			tot_flow += dinic_dfs(s, t, INF_FLOW);
		}
		return tot_flow;
	}

	std::vector<std::pair<flow_t, cost_t>> all_flows(int s, int t, cost_t max_cost = INF_COST - 1) {
		assert(s != t);
		build();
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
		build();
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
		int rev;
		flow_t cap;
	};
	struct input_edge_t {
		int from;
		int to;
		flow_t cap;
		flow_t rev_cap;
	};
	std::vector<input_edge_t> input;
	std::vector<int> adj_start;
	std::vector<edge_t> csr;
	std::vector<int> edge_pos;

	explicit Dinic(int N_) : N(N_) {}

	int add_edge(int from, int to, flow_t cap) {
		return add_bi_edge(from, to, cap, 0);
	}

	int add_bi_edge(int from, int to, flow_t cap, flow_t rev_cap) {
		assert(adj_start.empty());
		assert(cap >= 0);
		assert(rev_cap >= 0);
		input.push_back(input_edge_t{from, to, cap, rev_cap});
		return int(input.size()) - 1;
	}

	void build() {
		if (!adj_start.empty()) return;
		adj_start.assign(N+1, 0);
		for (const auto& e : input) {
			adj_start[e.from+1]++;
			adj_start[e.to+1]++;
		}
		for (int i = 0; i < N; i++) adj_start[i+1] += adj_start[i];
		csr.resize(2 * input.size());
		edge_pos.resize(input.size());
		std::vector<int> pos(adj_start.begin(), adj_start.end()-1);
		for (int i = 0; i < int(input.size()); i++) {
			const auto& e = input[i];
			int ef = pos[e.from]++, et = pos[e.to]++;
			csr[ef] = edge_t{e.to, et, e.cap};
			csr[et] = edge_t{e.from, ef, e.rev_cap};
			edge_pos[i] = ef;
		}
	}

	flow_t flow(int e) const {
		return input[e].cap - csr[edge_pos[e]].cap;
	}

	static constexpr tot_flow_t INF_FLOW = std::numeric_limits<tot_flow_t>::max() / 4;
	std::vector<int> buf;
	std::vector<int> level;
	std::vector<int> it;
	tot_flow_t dinic_dfs(int cur, int t, tot_flow_t f) {
		if (cur == t) return f;
		tot_flow_t cur_f = 0;
		assert(f > 0);
		for (; it[cur] < adj_start[cur+1]; it[cur]++) {
			edge_t& e = csr[it[cur]];
			int nxt = e.dest;
			if (level[nxt] == level[cur] + 1 && e.cap > 0) {
				flow_t v = flow_t(dinic_dfs(nxt, t, std::min<tot_flow_t>(f, e.cap)));
				e.cap -= v;
				csr[e.rev].cap += v;
				f -= v;
				cur_f += v;
				if (f == 0) break;
			}
		}
		return cur_f;
	}
	tot_flow_t dinic(int s, int t) {
		build();
		tot_flow_t tot_flow = 0;
		while (true) {
			buf.clear();
			buf.reserve(N);
			level.assign(N, -1);
			buf.push_back(s);
			level[s] = 0;
			for (int z = 0; z < int(buf.size()); z++) {
				int cur = buf[z];
				for (int a = adj_start[cur]; a < adj_start[cur+1]; a++) {
					const edge_t& e = csr[a];
					int nxt = e.dest;
					if (e.cap > 0 && level[nxt] == -1) {
						level[nxt] = level[cur] + 1;
						buf.push_back(nxt);
					}
				}
			}
			if (level[t] == -1) break;
			it.assign(adj_start.begin(), adj_start.end()-1);
			tot_flow += dinic_dfs(s, t, INF_FLOW);
		}
		return tot_flow;
	}
	tot_flow_t max_flow(int s, int t) { return dinic(s, t); }
};
