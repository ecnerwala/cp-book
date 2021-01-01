namespace pe::mod {

typedef long long ll;

ll invmod(ll a, ll m) {
	if (a == 1) return 1;
	return m - invmod(m % a, a) * m / a;
}

template<int MOD> struct Fp {
	static_assert(MOD >= 2);
	// static_assert(is_prime(MOD)); // we don't have this code

	int val;

	Fp() : val(0) {}
	Fp(int v) : val(v % MOD) { if (val < 0) val += MOD; }
	operator int () const { return val; }

	Fp inv() const {
		return Fp(int(invmod(val, MOD)));
	}

	Fp& operator += (const Fp &other) {
		val += other.val;
		if (val >= MOD) val -= MOD;
		return *this;
	}
	Fp& operator -= (const Fp &other) {
		val -= other.val;
		if (val < 0) val += MOD;
		return *this;
	}
	Fp& operator *= (const Fp &other) {
		val = int(ll(val) * other.val % MOD);
		return *this;
	}
	Fp& operator /= (const Fp &other) {
		assert(other.val != 0);
		return *this *= other.inv();
	}

	Fp operator + (const Fp &other) const {
		Fp res = *this;
		res += other;
		return res;
	}
	Fp operator - (const Fp &other) const {
		Fp res = *this;
		res -= other;
		return res;
	}
	Fp operator * (const Fp &other) const {
		Fp res = *this;
		res *= other;
		return res;
	}
	Fp operator / (const Fp &other) const {
		Fp res = *this;
		res /= other;
		return res;
	}

	Fp operator - () const {
		return Fp(MOD-val);
	}

	Fp operator + () const {
		return Fp(+val);
	}
};

template <typename T> T pow(T a, ll b) {
	T r = 1;
	while (b) {
		if (b & 1) r *= a;
		a *= a;
		b /= 2;
	}
	return r;
}

}
