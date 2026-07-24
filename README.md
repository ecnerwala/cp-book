# ecnerwala's CP Book

[![CI](https://github.com/ecnerwala/cp-book/actions/workflows/ci.yml/badge.svg)](https://github.com/ecnerwala/cp-book/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-github.io-blue?logo=github)](https://ecnerwala.github.io/cp-book/)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fecnerwala.github.io%2Fcp-book%2Fcoverage%2Fbadge.json)](https://ecnerwala.github.io/cp-book/coverage/)

This is my library of reference code for competitive programming. The goal is to
write generic, fast, and clean algorithm implementations for use in contests
like CodeForces or ICPC.

## Building

Build using

```sh
cmake -B build
cmake --build build
```

Test with

```sh
ctest --test-dir build
```

or directly with

```sh
./build/tests
```

## Library Checker verification

Solutions to [Library Checker](https://judge.yosupo.jp/) problems live in
`verify/` and are verified in CI with
[competitive-verifier](https://github.com/competitive-verifier/competitive-verifier);
results and browsable (bundled) code are hosted at
https://ecnerwala.github.io/cp-book/.

The verify files include library headers directly (e.g. `#include "fft.hpp"`,
resolved from `src/`). To produce a single submittable file with all headers
inlined:

```sh
scripts/bundle.py verify/convolution_mod.test.cpp > submission.cpp
```

(uses [uv](https://docs.astral.sh/uv/); alternatively run with python3 after
`pip install competitive-verifier`.)

To run verification locally:

```sh
uvx competitive-verifier oj-resolve --include src verify --exclude third_party \
    --config .competitive-verifier/config.toml > verify_files.json
uvx competitive-verifier verify --verify-json verify_files.json
```

## License and Attribution

All code in this book is written by me and CC0 licensed unless otherwise noted
in the file. Inspiration is largely drawn from KACTL
(https://github.com/kth-competitive-programming/kactl/) and other references.
