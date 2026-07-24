# ecnerwala's Competitive Programming Book

[![CI](https://github.com/ecnerwala/cp-book/actions/workflows/ci.yml/badge.svg)](https://github.com/ecnerwala/cp-book/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-github.io-blue?logo=github)](https://ecnerwala.github.io/cp-book/)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fecnerwala.github.io%2Fcp-book%2Fcoverage%2Fbadge.json)](https://ecnerwala.github.io/cp-book/coverage/)

This is my library of reference code for competitive programming. The goal is to
write generic, fast, and clean algorithm implementations for use in contests
like CodeForces or ICPC.

The library lives in `src/` as standalone headers (`#include "fft.hpp"`,
compiled with `-I src`). Browsable source, verification results, and coverage
are hosted at https://ecnerwala.github.io/cp-book/.

## Building and testing

```sh
cmake -B build && cmake --build build
ctest --test-dir build          # or ./build/tests
```

## Library Checker verification

`verify/` holds [Library Checker](https://judge.yosupo.jp/) solutions,
verified in CI with
[competitive-verifier](https://github.com/competitive-verifier/competitive-verifier).

```sh
# Single submittable file with headers inlined
scripts/bundle.py verify/convolution_mod.test.cpp > submission.cpp

# Run verification locally
uvx competitive-verifier oj-resolve --include src verify --exclude third_party \
    --config .competitive-verifier/config.toml > verify_files.json
uvx competitive-verifier verify --verify-json verify_files.json
```

## License and attribution

All code is written by me and CC0 licensed unless otherwise noted in the
file. Inspiration is largely drawn from
[KACTL](https://github.com/kth-competitive-programming/kactl/) and other
references.
