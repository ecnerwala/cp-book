#pragma once
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,popcnt,abm,mmx") // Safe for yandex

#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,bmi,bmi2,mmx,avx,avx2,fma") // Requires AVX2

// See https://codeforces.com/blog/entry/96344

inline void disable_denormal_floats() {
	// https://stackoverflow.com/a/8217313
	#define CSR_FLUSH_TO_ZERO         (1 << 15)
	unsigned csr = __builtin_ia32_stmxcsr();
	csr |= CSR_FLUSH_TO_ZERO;
	__builtin_ia32_ldmxcsr(csr);
	#undef CSR_FLUSH_TO_ZERO
}
