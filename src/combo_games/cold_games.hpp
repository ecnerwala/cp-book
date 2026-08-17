#pragma once

#include <iostream>
#include <vector>
#include <optional>
#include <array>
#include <cassert>

namespace wala {


// TODO: Make this generic over numerator type, e.g. bignum
struct dyadic {
	int n = 0, d = 0;
	friend std::ostream& operator << (std::ostream& o, dyadic v) {
		return o << v.n << "/2^" << v.d;
	}
	friend auto operator == (dyadic a, dyadic b) {
		return a.n == b.n && a.d == b.d;
	}
	friend auto operator <=> (dyadic a, dyadic b) {
		return a.n << std::max(b.d - a.d, 0) <=> b.n << std::max(a.d - b.d, 0);
	}
	dyadic operator - () const {
		return {-n, d};
	}
	friend dyadic operator + (dyadic a, dyadic b) {
		int nn = (a.n << std::max(b.d - a.d, 0)) + (b.n << std::max(a.d - b.d, 0));
		int nd = std::max(a.d, b.d);
		if (!nn) nd = 0;
		else {
			int x = std::min(__builtin_ctz(nn), nd);
			nn >>= x;
			nd -= x;
		}
		return {nn, nd};
	}
};

// TODO: Make this generic over dyadic type
struct cold_ish {
	dyadic v;
	bool star = false;
	friend std::ostream& operator << (std::ostream& o, cold_ish d) {
		return o << d.v << (d.star ? "*" : "");
	}

	friend bool operator == (cold_ish a, cold_ish b) = default;

	cold_ish operator - () const {
		return {-v, star};
	}
	friend cold_ish operator + (cold_ish a, cold_ish b) {
		return {a.v + b.v, bool(a.star ^ b.star)};
	}
};

// If L_options <| a <| R_options and no (direct) ancestors of a
// satisfy this (by breaking it on the correct side), then G = a.
//
// Proof:
// L can move to L_options - a <| 0 which means R wins as the next player,
// or G - a.R so R moves to some R_options - a.R <= 0 which means R wins as the 2nd player
// (this requires that a.R breaks it by being >= some R_options and not <= some L_options,
// which holds if everything here is cold-ish).

// 0 are moves for l
// TODO: Could optimize
inline std::optional<cold_ish> cold_ish_game(std::array<std::vector<cold_ish>, 2> moves) {
	std::array<std::optional<cold_ish>, 2> best;
	for (int z = 0; z < 2; z++) {
		std::optional<cold_ish> v;
		for (auto m : moves[z]) {
			if (z) m = -m;
			// no star is "bigger" in this sense
			if (!v || m.v > v->v || (m.v == v->v && m.star < v->star)) {
				v = m;
			}
		}
		if (z && v) v = -*v;
		best[z] = v;
	}
	auto fuzzy_less = [&](std::optional<cold_ish> a, std::optional<cold_ish> b) -> bool {
		if (!a) return true;
		if (!b) return true;
		return (a->v < b->v) || (a->v == b->v && a->star != b->star);
	};
	if (best[0] && best[1] && fuzzy_less(*best[1], *best[0])) {
		// invalid
		return std::nullopt;
	}

	// best[0] <= best[1], so we should be able to find someone here

	// the only time we return star is if the lower bound == the upper bound
	if (best[0] && best[1] && *best[0] == *best[1]) {
		auto cnd = *best[0];
		cnd.star ^= true;
		assert(fuzzy_less(best[0], cnd) && fuzzy_less(cnd, best[1]));
		return cnd;
	}
	assert(!best[0] || !best[1] || best[0]->v < best[1]->v);

	cold_ish zero{{0}, false};
	bool le_zero = fuzzy_less(best[0], zero);
	bool ge_zero = fuzzy_less(zero, best[1]);
	assert(le_zero || ge_zero);
	if (le_zero && ge_zero) {
		return zero;
	}
	bool flip = le_zero;
	if (flip) {
		std::swap(le_zero, ge_zero);
		std::swap(best[0], best[1]);
		if (best[0]) best[0] = -*best[0];
		if (best[1]) best[1] = -*best[1];
	}
	assert(ge_zero);
	assert(!le_zero);
	// check if there's an integer that's good
	assert(best[0]);

	// best[0] >= 0
	int int_cnd = (best[0]->star ? best[0]->v.n - 1 : best[0]->v.n);
	assert(int_cnd >= 0);
	int_cnd >>= best[0]->v.d;
	int_cnd++;
	cold_ish cnd = {dyadic(int_cnd)};
	assert(fuzzy_less(best[0], cnd));
	if (fuzzy_less(cnd, best[1])) {
		return flip ? -cnd : cnd;
	}

	// otherwise, cnd is too big, and we're just doing the dyadic rational thing
	assert(best[1]);
	// extra 1 to allow the middle
	int d = std::max(best[0]->v.d, best[1]->v.d)+1;
	// exclusive bounds
	int v0 = (best[0]->v.n << (d - best[0]->v.d)) - best[0]->star;
	int v1 = (best[1]->v.n << (d - best[1]->v.d)) + best[1]->star;
	assert(v1 - v0 >= 2);
	int b = 31 - __builtin_clz((v1-1) ^ v0);
	cnd = cold_ish{dyadic((v1-1) >> b, d-b), false};
	assert(cnd.v.n & 1);
	assert(fuzzy_less(best[0], cnd) && fuzzy_less(cnd, best[1]));
	return flip ? -cnd : cnd;
}

} // namespace wala
