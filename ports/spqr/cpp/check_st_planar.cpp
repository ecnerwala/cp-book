// Empirical check of the st-order and planarity properties of planar_spqr_tree, to pin down the Lean spec.
#include "graph/spqr_tree.hpp"
#include <iostream>
#include <map>
#include <set>
#include <functional>
#include <numeric>
using namespace std;
using T = wala::planar_spqr_tree; using NT = T::node_type;
int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { if (fails < 20) cerr << "FAIL " << msg << '\n'; fails++; } } while (0)
int main() {
	int NV, NE, tern; cin >> NV >> NE >> tern;
	vector<array<int,2>> edges(NE); for (auto& [u,v] : edges) cin >> u >> v;
	int K; cin >> K; vector<int> vo(K); for (auto& x : vo) cin >> x;
	int L; cin >> L; vector<int> eo(L); for (auto& x : eo) cin >> x;
	auto t = T::build(NV, edges, tern, vo, eo);
	int n = t.size();
	int cntR = 0, cntNonplanar = 0, cntDense = 0;
	for (int i = 0; i < n; i++) {
		auto ty = t.types[i];
		int s = t.node_nvs.bounds[i], e = t.node_nvs.bounds[i+1];
		int es = t.node_nes.bounds[i], ee = t.node_nes.bounds[i+1];
		// (1) st-numbering: every interior nv has a lower and a higher neighbour; edges oriented low->high
		if (ty == NT::S || ty == NT::R || ty == NT::P) {
			vector<bool> lo(e - s), hi(e - s);
			for (int ne = es; ne < ee; ne++) {
				auto [a, b] = t.node_edges[ne].nvs;
				CHECK(a < b, "edge orientation node " << i);
				hi[a - s] = true; lo[b - s] = true;
			}
			for (int nv = s + 1; nv + 1 < e; nv++) CHECK(lo[nv - s] && hi[nv - s], "st-numbering node " << i << " nv " << nv);
		}
		// (2) dominance order among non-cap edges (and whether the cap breaks it)
		bool hascap = ty != NT::F && ty != NT::V && !(ty == NT::Q && t.ch.bounds[i+1] > t.ch.bounds[i]);
		for (int a = es + hascap; a < ee; a++) for (int b = a + 1; b < ee; b++) {
			auto pa = t.node_edges[a].nvs, pb = t.node_edges[b].nvs;
			if (pa == pb) continue;
			CHECK(!(pb[0] <= pa[0] && pb[1] <= pa[1]), "dominance node " << i << " " << a << " " << b);
		}
		// (3) bracket adjacency rows
		for (int nv = s; nv < e; nv++) {
			for (int side = 0; side < 2; side++) {
				int lo_ = t.node_adj.bounds[2*nv+side], hi_ = t.node_adj.bounds[2*nv+side+1];
				for (int k = lo_; k < hi_; k++) {
					int d = t.node_adj.dat[k].dest_nv;
					if (e - s > 1) CHECK(side == 0 ? d < nv : d > nv, "adj side node " << i << " nv " << nv);
					if (k + 1 < hi_) CHECK(t.node_adj.dat[k+1].dest_nv <= d, "adj order node " << i << " nv " << nv << " side " << side);
				}
			}
		}
		// (4) per-node embedding
		if (ty == NT::S || ty == NT::R || ty == NT::P) {
			if (ty == NT::R) cntR++;
			if (!t.node_planar[i]) {
				cntNonplanar++;
				CHECK(ty == NT::R, "nonplanar non-R node " << i);
				int V = e - s, E = ee - es;
				if (E > 3 * V - 6) cntDense++;
				cout << "NP " << V; for (int ne = es; ne < ee; ne++) cout << ' ' << t.node_edges[ne].nvs[0]-s << ' ' << t.node_edges[ne].nvs[1]-s; cout << '\n';
				for (int q = 4*es; q < 4*ee; q++) CHECK(t.ne_embedding.rot_adj[q] == -1, "nonplanar rot_adj set");
				continue;
			}
			auto& ra = t.ne_embedding.rot_adj;
			auto vert_of = [&](int q) { return t.node_edges[q >> 2].nvs[(q >> 1) & 1]; };
			for (int q = 4*es; q < 4*ee; q++) {
				int r = ra[q];
				CHECK(r >= 4*es && r < 4*ee, "rot_adj range node " << i);
				CHECK(ra[r] == q, "rot_adj involution node " << i);
				CHECK((r & 1) != (q & 1), "rot_adj dir bit node " << i);
				CHECK(vert_of(r) == vert_of(q), "rot_adj same vertex node " << i);
			}
			// vertex orbits: q -> ra[q ^ 1] must visit exactly the quarter-edges at that vertex (one cycle per vertex)
			{
				vector<bool> seen(4*(ee-es));
				int cycles = 0;
				for (int q = 4*es; q < 4*ee; q++) if (!seen[q-4*es]) {
					cycles++;
					int x = q; do { seen[x-4*es] = true; x = ra[x ^ 1]; } while (x != q);
				}
				CHECK(cycles == 2 * (e - s), "vertex rotation cycles node " << i << " got " << cycles << " want " << 2*(e-s));
			}
			// faces: q -> ra[q ^ 3]
			{
				vector<bool> seen(4*(ee-es));
				int F = 0;
				for (int q = 4*es; q < 4*ee; q++) if (!seen[q-4*es]) {
					F++;
					int x = q; do { seen[x-4*es] = true; x = ra[x ^ 3]; } while (x != q);
				}
				int V = e - s, E = ee - es;
				CHECK(F == 2 * (2 - V + E), "euler node " << i << " V=" << V << " E=" << E << " F/2=" << F/2.0);
			}
		}
	}
	// (5) glued embedding
	auto emb = wala::planar_embed(t);
	bool allplanar = all_of(t.node_planar.begin(), t.node_planar.end(), identity{});
	CHECK(emb.has_value() == allplanar, "planar_embed presence");
	if (emb) {
		auto& ra = emb->rot_adj;
		auto vert_of = [&](int q) { return edges[q >> 2][(q >> 1) & 1]; };
		for (int q = 0; q < 4*NE; q++) {
			int r = ra[q];
			CHECK(r >= 0 && r < 4*NE, "glued range");
			if (r < 0) continue;
			CHECK(ra[r] == q, "glued involution");
			CHECK((r & 1) != (q & 1), "glued dir bit");
			CHECK(vert_of(r) == vert_of(q), "glued same vertex");
		}
		vector<int> deg(NV); for (auto [u,v] : edges) { deg[u]++; deg[v]++; }
		vector<bool> seen(4*NE); int cyc = 0;
		for (int q = 0; q < 4*NE; q++) if (!seen[q]) { cyc++; int x = q; do { seen[x] = true; x = ra[x ^ 1]; } while (x != q); }
		int nonisolated = 0; for (int v = 0; v < NV; v++) nonisolated += deg[v] > 0;
		CHECK(cyc == 2 * nonisolated, "glued vertex cycles " << cyc << " vs " << 2*nonisolated);
		fill(seen.begin(), seen.end(), false); int F = 0;
		for (int q = 0; q < 4*NE; q++) if (!seen[q]) { F++; int x = q; do { seen[x] = true; x = ra[x ^ 3]; } while (x != q); }
		// components (ignoring isolated vertices)
		vector<int> p(NV); iota(p.begin(), p.end(), 0);
		function<int(int)> f = [&](int x) { return p[x] == x ? x : p[x] = f(p[x]); };
		for (auto [u,v] : edges) p[f(u)] = f(v);
		int C = 0; for (int v = 0; v < NV; v++) C += deg[v] > 0 && f(v) == v;
		// F counts each face twice (both dirs); each component embedded separately: V - E + F = 2C
		CHECK(F == 2 * (2*C - nonisolated + NE), "glued euler F/2=" << F/2.0 << " V=" << nonisolated << " E=" << NE << " C=" << C);
	}
	cout << "items=" << n << " R=" << cntR << " nonplanarR=" << cntNonplanar << " (dense " << cntDense << ") glued=" << int(emb.has_value()) << " fails=" << fails << '\n';
	return fails ? 1 : 0;
}
