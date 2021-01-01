#include <bits/stdc++.h>
#include "fft.hpp"
#include "modnum.hpp"
using namespace std;

using ntt1_t = fft<21, modnum<998244353>>;
ntt1_t ntt1(421);
ntt1_t::arr a;
ntt1_t::arr b;

vector<int> multiply_blocked(const vector<int>& s, const vector<int>& t) {
	vector<long long> res(s.size() + t.size());
	const int B = (1 << 20);
	for (int i = 0; i < int(s.size()); i += B) {
		memset(a, 0, sizeof(a));
		for (int x = 0; x < B && i + x < int(s.size()); x++) {
			a[x] = s[i + x];
		}
		ntt1(a);
		for (int j = 0; j < int(t.size()); j += B) {
			memset(b, 0, sizeof(b));
			for (int y = 0; y < B && j + y < int(t.size()); y++) {
				b[y] = t[j + y];
			}
			ntt1(b);
			for (int x = 0; x < ntt1.N; x++) {
				b[x] *= a[x];
			}
			ntt1(b, true);
			for (int x = 0; x < ntt1.N; x++) {
				if (x + i + j < int(res.size())) {
					res[x + i + j] += int(b[x]);
				} else {
					assert(b[x] == 0);
				}
			}
		}
	}
	vector<int> r;
	long long carry = 0;
	for (int i = 0; i < int(res.size()); i++) {
		carry += res[i];
		r.push_back(int(carry % 10));
		carry /= 10;
	}
	assert(carry == 0);
	while (r.back() == 0) {
		r.pop_back();
	}
	return r;
}

using ntt2_t = fft<24, modnum<754974721>>;
ntt2_t ntt2(362);
ntt2_t::arr a1;
ntt2_t::arr b1;
using ntt3_t = fft<24, modnum<469762049>>;
ntt3_t ntt3(320192759);
ntt3_t::arr a2;
ntt3_t::arr b2;

vector<int> multiply(const vector<int>& s, const vector<int>& t) {
	assert(s.size() + t.size() <= ntt2.N);
	assert(s.size() + t.size() <= ntt3.N);
	memset(a1, 0, sizeof(a1));
	memset(a2, 0, sizeof(a2));
	memset(b1, 0, sizeof(b1));
	memset(b2, 0, sizeof(b2));
	for (int i = 0; i < int(s.size()); i++) {
		a1[i] = s[i];
		a2[i] = s[i];
	}
	for (int i = 0; i < int(t.size()); i++) {
		b1[i] = t[i];
		b2[i] = t[i];
	}
	ntt2(a1);
	ntt2(b1);
	for (int i = 0; i < ntt2.N; i++) {
		b1[i] *= a1[i];
	}
	ntt2(b1, true);

	ntt3(a2);
	ntt3(b2);
	for (int i = 0; i < ntt3.N; i++) {
		b2[i] *= a2[i];
	}
	ntt3(b2, true);

	vector<int> res;
	long long carry = 0;
	long long m1 = ntt2_t::num_t::MOD;
	for (int i = 0; i < int(s.size() + t.size()); i++) {
		long long v1 = int(b1[i]);
		long long v2 = int(b2[i]);
		ntt3_t::num_t num(v2 - v1);
		num /= ntt3_t::num_t(m1);
		long long v = v1 + int(num) * m1;
		carry += v;
		res.push_back(carry % 10);
		carry /= 10;
	}
	assert(carry == 0);
	while (res.back() == 0) {
		res.pop_back();
	}
	return res;
}

using cfft_t = fft<24, complex<double>>;
cfft_t cfft(std::polar(1., 2 * M_PI / (1 << 24)));
cfft_t::arr ac;
cfft_t::arr bc;

vector<int> multiply_complex(const vector<int>& s, const vector<int>& t) {
	assert(s.size() + t.size() <= cfft.N);
	memset(ac, 0, sizeof(ac));
	memset(bc, 0, sizeof(bc));
	for (int i = 0; i < int(s.size()); i++) {
		ac[i] = s[i];
	}
	for (int i = 0; i < int(t.size()); i++) {
		bc[i] = t[i];
	}
	cerr << "HI\n";
	cfft(ac);
	cerr << "HI\n";
	cfft(bc);
	cerr << "HI\n";
	for (int i = 0; i < ntt2.N; i++) {
		bc[i] *= ac[i];
	}
	cerr << "HI\n";
	cfft(bc, true);
	cerr << "HI\n";

	vector<int> res;
	long long carry = 0;
	for (int i = 0; i < int(s.size() + t.size()); i++) {
		long long v = (long long)(bc[i].real() + 0.5);
		carry += v;
		res.push_back(carry % 10);
		carry /= 10;
	}
	cerr << "HI\n";
	assert(carry == 0);
	while (!res.empty() && res.back() == 0) {
		res.pop_back();
	}
	cerr << "HI\n";
	return res;
}

int main() {
	ios_base::sync_with_stdio(0);
	/*
	for (int i = 0; i < ntt1.N / 2; i++) {
		a[i] = i;
	}
	ntt1(a);
	ntt1(a, true);
	for (int i = 0; i < ntt1.N / 2; i++) {
		b[i] = 1;
	}
	ntt1(a);
	ntt1(b);
	for (int i = 0; i < ntt1.N; i++) {
		a[i] *= b[i];
	}
	ntt1(a, true);
	for (int i = 0; i < 10; i++) {
		cout << a[i] << ' ';
	}
	cout << '\n';
	*/

	int Q; cin >> Q;
	for (int q = 0; q < 7; q++) {
		cerr << "Case #" << q << ":" << '\n';
		string s, t; cin >> s >> t;
		cerr << s.size() << ' ' << t.size() << '\n';
		vector<int> sn, tn;
		for (char c : s) {
			sn.push_back(c - '0');
		}
		for (char c : t) {
			tn.push_back(c - '0');
		}
		reverse(sn.begin(), sn.end());
		reverse(tn.begin(), tn.end());
		//vector<int> res = multiply_blocked(sn, tn);
		//vector<int> res = multiply(sn, tn);
		vector<int> res = multiply_complex(sn, tn);
		string r;
		for (int i : res) {
			r.push_back(char(i + '0'));
		}
		reverse(r.begin(), r.end());
		cout << r << '\n';
	}

}
