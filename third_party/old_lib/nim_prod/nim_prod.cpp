#include<bits/stdc++.h>
using namespace std;

using ll = unsigned long long;

/**
 * Author: Andrew He
 * Description: Nim Product.
 */
ll _nimProd2[64][64];
ll nimProd2(int i, int j) {
  if (_nimProd2[i][j]) return _nimProd2[i][j];
  if ((i & j) == 0) return _nimProd2[i][j] = 1ull << (i|j);
  int a = (i&j) & -(i&j);
  return _nimProd2[i][j] = nimProd2(i ^ a, j) ^ nimProd2((i ^ a) | (a-1), (j ^ a) | ((i|j) & (a-1)));
}
ll nimProd(ll x, ll y) {
  ll res = 0;
	for (int i = 0; x >> i; i++)
    if ((x >> i) & 1)
      for (int j = 0; y >> j; j++)
        if ((y >> j) & 1)
          res ^= nimProd2(i, j);
  return res;
}

int main() {
	ios_base::sync_with_stdio(0), cin.tie(0), cout.tie(0);
	for (int i = 0; i < 64; i++) {
		for (int j = 0; j < 64; j++) {
			ll v = nimProd2(i, j);
			cout << i << ' ' << j << ' ' << v << '\n';
		}
	}

	// inverses?
	for (int i = 1; i < 64; i++) {
	  for (int j = 1; j < 64; j++) {
	    if (nimProd(i, j) == 1) {
	      cerr << i << ' ' << j << '\n';
	    }
	  }
	}

	return 0;
}
