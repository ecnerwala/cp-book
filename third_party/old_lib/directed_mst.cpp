#include<bits/stdc++.h>
using namespace std;

mt19937 mt(48);

const int MAXN = 1.1e6;
const int MAXM = 1.1e6;

int N, M;
int A[MAXN];

struct node {
	node* ch[2];
	int val;
	int target;
};
const int MAXNODES = MAXN + 2 * MAXM;
int NODES = 0;
node node_pool[MAXNODES];

node* edges[MAXN];

node* meld(node* a, node* b) {
	if (!a) return b;
	if (!b) return a;
	if (b->val > a->val) {
		swap(a, b);
	}
	b->val -= a->val;
	int x = mt() & 1;
	a->ch[x] = meld(b, a->ch[x]);
	return a;
}

int par[MAXN];
int getpar(int a) { return par[a] == -1 ? a : (par[a] = getpar(par[a])); }

bool inStack[MAXN];

int main() {
	ios_base::sync_with_stdio(0), cin.tie(0), cout.tie(0);
	cin >> N >> M;
	for (int i = 1; i <= N; i++) cin >> A[i];

	auto addEdge = [](int source, int target, int v) {
		node* n = &node_pool[NODES++];
		n->val = v;
		n->ch[0] = n->ch[1] = nullptr;
		n->target = target;
		edges[source] = meld(edges[source], n);
	};

	for (int i = 1; i <= N; i++) {
		addEdge(i, 0, 0);
	}
	for (int e = 0; e < M; e++) {
		int x, y; cin >> x >> y;
		addEdge(x, y, A[y]);
		addEdge(y, x, A[x]);
	}

	memset(par, -1, sizeof(par));

	long long ans = 0;
	for (int i = 1; i <= N; i++) {
		if (getpar(i) != i) continue;

		stack<int> st;
		st.emplace(0);
		inStack[0] = true;

		int cur = i;
		while (!st.empty()) {
			st.push(cur);
			inStack[cur] = true;
			int nxt = getpar(edges[cur]->target);
			ans += edges[cur]->val;
			edges[cur] = meld(edges[cur]->ch[0], edges[cur]->ch[1]);
			cur = nxt;
			while (inStack[cur]) {
				int o = st.top(); st.pop();
				inStack[o] = false;
				if (o != cur) {
					par[o] = cur;
					edges[cur] = meld(edges[o], edges[cur]);
				}
			}
		}
	}
	cout << ans << '\n';

	return 0;
}
