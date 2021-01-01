#include <bits/stdc++.h>
#include "nt/arith_sum.hpp"
#include "modnum.hpp"
using namespace std;

using namespace nt::arith_sum;

using ll = long long;
using num = modnum<ll(1e9)+7>;


div_vector<ll(1e6), num> dv, prod, divi, rt;

int main() {
	ios_base::sync_with_stdio(0);

	dv.fill([](ll v) -> num { return num(v) * num(v+1) / 2; });
	prod = convolve(dv, dv);
	cerr << "Product done\n";

	divi = divide(prod, dv);
	for (auto it : divi) {
		if (it.second != dv[it.first]) {
			cerr << it.first << ' ' << it.second << ' ' << dv[it.first] << '\n';
			assert(it.second == dv[it.first]);
		}
	}

	rt = sqrt(prod);
	for (auto it : rt) {
		if (it.second != dv[it.first]) {
			cerr << it.first << ' ' << it.second << ' ' << dv[it.first] << '\n';
			assert(it.second == dv[it.first]);
		}
	}

	auto tmp = div_vector<ll(1e6), num>::POLY(1);
	for (auto it : tmp) {
		if (it.second != dv[it.first]) {
			cerr << it.first << ' ' << it.second << ' ' << dv[it.first] << '\n';
			assert(it.second == dv[it.first]);
		}
	}

	cout << "TEST SUCCEEDED\n";
	return 0;
}
