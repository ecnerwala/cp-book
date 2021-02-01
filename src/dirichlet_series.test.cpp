#include "dirichlet_series.hpp"

#include "catch.hpp"
#include <bits/stdc++.h>

#include "modnum.hpp"

namespace dirichlet_series {

namespace test {

div_vector_layout layout;
template <typename T> using dv_values = dirichlet_series_values<layout, T>;
template <typename T> using dv_prefix = dirichlet_series_prefix<layout, T>;

template <typename T>
dv_values<T> multiply_slow(const dv_values<T>& a, const dv_values<T>& b) {
	dv_values<T> r;
	for (int i = 1; i < layout.len; i++) {
		for (int j = 1; j < layout.len; j++) {
			int k = layout.get_value_bucket(layout.get_bucket_bound(i) * layout.get_bucket_bound(j));
			if (k < layout.len) r.st[k] += a.st[i] * b.st[j];
		}
	}
	return r;
}

TEMPLATE_TEST_CASE("Dirichlet series multiplication and inverse", "[dirichlet]", modnum<int(1e9)+7>, int64_t) {
	using num = TestType;
	for (int N = 1; N <= 30; N++) {
		INFO("N = " << N);
		std::mt19937 mt(48);
		layout = div_vector_layout(N);
		dv_prefix<num> a([&](int64_t x) { return num(x); });
		dv_prefix<num> b([&](int64_t x) { return num(x) * num(x+1) / num(2); });
		dv_prefix<num> slow_res(multiply_slow(dv_values<num>(a), dv_values<num>(b)));
		dv_prefix<num> fast_res = a * b;
		for (int i = 1; i < layout.len; i++) {
			INFO("i = " << i);
			REQUIRE(slow_res.st[i] == fast_res.st[i]);
		}
		dv_prefix<num> a_2 = fast_res / b;
		for (int i = 1; i < layout.len; i++) {
			INFO("i = " << i);
			REQUIRE(a.st[i] == a_2.st[i]);
		}
		dv_prefix<num> b_2 = fast_res / a;
		for (int i = 1; i < layout.len; i++) {
			INFO("i = " << i);
			REQUIRE(b.st[i] == b_2.st[i]);
		}
		if constexpr (!std::is_same_v<num, int64_t>) {
			dv_prefix<num> rt_a = sqrt(a);
			dv_prefix<num> a_3 = rt_a * rt_a;
			for (int i = 1; i < layout.len; i++) {
				INFO("i = " << i);
				REQUIRE(a.st[i] == a_3.st[i]);
			}
			dv_prefix<num> rt_b = sqrt(b);
			dv_prefix<num> b_3 = rt_b * rt_b;
			for (int i = 1; i < layout.len; i++) {
				INFO("i = " << i);
				REQUIRE(b.st[i] == b_3.st[i]);
			}
		}
	}
}

TEMPLATE_TEST_CASE("Dirichlet series euler transform", "[dirichlet]", modnum<int(1e9)+7>) {
	using num = TestType;
	for (int N : {1, 2, 3, 4, 5, 24, 25, 26, 99, 100, 101}) {
		INFO("N = " << N);
		std::mt19937 mt(48);
		layout = div_vector_layout(N);
		dv_prefix<num> a([&](int64_t x) { return num(x); });
		dv_values<num> primes = inverse_euler_transform(a);
		dv_values<num> primes_slow;
		for (int v = 2; v <= N; v++) {
			bool is_prime = true;
			for (int p = 2; p * p <= v; p++) {
				if (v % p == 0) {
					is_prime = false;
					break;
				}
			}
			primes_slow[v] += is_prime;
		}
		for (int i = 1; i < layout.len; i++) {
			REQUIRE(primes_slow[i] == primes[i]);
		}
	}
}

}

}
