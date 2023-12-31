#pragma once

#include <vector>
#include <cassert>

#include "yc.hpp"

// Direct a graph into a DAG so that given source and sink are the unique sources/sinks.
// If there are any biconnected components not on the path from the source to
// the sink, they will not be output, modify the code if necessary.
// Returns a topological sort. Back out the edge directions yourself.
inline std::vector<int> make_st_dag(const std::vector<std::vector<int>>& adj, int source = -1, int sink = -1) {
	int N = int(adj.size());

	// What's even going on lol
	if (N == 0) return {};

	// Make some arbitrary choices as defaults
	if (source == -1 && sink == -1) source = 0;
	if (source == -1) source = adj[sink].empty() ? sink : adj[sink][0];
	if (sink == -1) sink = adj[source].empty() ? source : adj[source][0];

	std::vector<int> depth(N, -1);
	std::vector<int> lowval(N);
	std::vector<bool> has_sink(N);
	std::vector<std::vector<int>> ch(N);
	std::y_combinator([&](auto self, int cur, int prv) -> void {
		depth[cur] = prv != -1 ? depth[prv] + 1 : 0;
		lowval[cur] = depth[cur];
		ch[cur].reserve(adj[cur].size());
		has_sink[cur] = (cur == sink);
		for (int nxt : adj[cur]) {
			if (nxt == prv) continue;
			if (depth[nxt] == -1) {
				ch[cur].push_back(nxt);
				self(nxt, cur);
				lowval[cur] = std::min(lowval[cur], lowval[nxt]);
				if (has_sink[nxt]) has_sink[cur] = true;
			} else if (depth[nxt] < depth[cur]) {
				lowval[cur] = std::min(lowval[cur], depth[nxt]);
			} else {
				// down edge
			}
		}
	})(source, -1);

	// true is after, false is before
	std::vector<bool> edge_dir(N, false);
	std::vector<int> lst_nxt(N, -1);
	auto lst = std::y_combinator([&](auto self, int cur) -> std::array<int, 2> {
		std::array<int, 2> res{cur, cur};
		for (int nxt : ch[cur]) {
			// If we're on the path to the sink, mark it as downwards.

			// Comment out this line to direct extra bcc's as extra sinks
			if (!has_sink[nxt] && lowval[nxt] >= depth[cur]) continue;

			bool d = (has_sink[nxt] || lowval[nxt] >= depth[cur]) ? true : !edge_dir[lowval[nxt]];
			edge_dir[depth[cur]] = d;

			auto ch_res = self(nxt);

			// Join res and ch
			if (!d) std::swap(res, ch_res);
			lst_nxt[std::exchange(res[1], ch_res[1])] = ch_res[0];
		}
		return res;
	})(source);

	std::vector<int> res; res.reserve(N);
	int cur = lst[0];
	while (cur != -1) {
		res.push_back(cur);
		cur = lst_nxt[cur];
	}
	return res;
}
