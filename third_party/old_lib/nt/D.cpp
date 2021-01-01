#include <bits/stdc++.h>
using namespace std;

using ll = long long;

constexpr ll floor_sqrt(ll n) {
	if (n == 0) return 0; assert(n > 0);
	return std::llround(std::trunc(std::sqrt(n)));
}

const int X = 1e8;

bitset<X> is_prime;
vector<int> pr;

int mu[X];
ll small_mu_sum[X];

void init(){
	is_prime.flip();
	is_prime[0] = is_prime[1] = false;
	mu[1] = 1;
	for(int i = 2; i < X; i++){
		if(is_prime[i]){
			pr.push_back(i);
			mu[i] = -1;
		}
		for(int p : pr){
			if(ll(i) * p >= X) break;
			is_prime[i * p] = false;
			mu[i * p] = -mu[i];
			if(i % p == 0){
				mu[i * p] = 0;
				break;
			}
		}
	}
	small_mu_sum[0] = 0;
	for(int i = 1; i < X; i++) small_mu_sum[i] = small_mu_sum[i-1] + mu[i];
}

void solve(){
	ll N;
	cin >> N;
	ll ans = 0;

	vector<bool> done(N / X + 1, false);
	vector<ll> memo(N / X + 1, 0);

	function<ll(ll)> mu_sum;
	mu_sum = [&](ll n) -> ll {
		if(n < X) return small_mu_sum[n];
		ll idx = N / n;
		if(!done[idx]){
			ll res = 1;
			ll f = 2;
			while(f <= n){
				ll r = n / f;
				ll g = n / r;
				res -= (g - f + 1) * mu_sum(r);
				f = g + 1;
			}
			memo[idx] = res;
			done[idx] = true;
		}
		return memo[idx];
	};

	ll f = 1;
	while(f <= N){
		ll g = N / (N / f);
		ans += (mu_sum(g) - mu_sum(f - 1)) * floor_sqrt(N / g);
		f = g + 1;
	}
	cout << ans << '\n';
}

int main(){
	ios_base::sync_with_stdio(0), cin.tie(0), cout.tie(0);
	init();
	solve();
}