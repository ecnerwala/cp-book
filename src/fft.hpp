#pragma once

#include "fft/core.hpp"
#include "fft/engines/ntt.hpp"
#include "fft/engines/real.hpp"
#include "fft/engines/split.hpp"
#include "fft/engines/crt.hpp"
#include "fft/engines/wrappers.hpp"
#include "fft/multiply.hpp"
#include "fft/series.hpp"
#include "fft/poly.hpp"
#include "fft/online.hpp"

namespace ecnerwala {
namespace fft {
namespace engines {

static_assert(engine<matrix<ntt<modnum<998244353>>, 2>>);
static_assert(engine<trunc<ntt<modnum<998244353>>, 3>>);
// tracked inner engines work when the accumulated scale fits the budget (N <= 2)
static_assert(engine<matrix<split<modnum<int(1e9)+7>>, 2>>);
static_assert(engine<trunc<crt<modnum<int(1e9)+7>>, 2>>);
// the stable variants keep tracked inner engines sound at any N
static_assert(engine<matrix_stable<split<modnum<int(1e9)+7>>, 3>>);
static_assert(engine<trunc_stable<crt<modnum<int(1e9)+7>>, 3>>);

/* namespace engines */ }
/* namespace fft */ }
/* namespace ecnerwala */ }
