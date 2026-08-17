# ecnerwala's Competitive Programming Book

[![CI](https://github.com/ecnerwala/cp-book/actions/workflows/ci.yml/badge.svg)](https://github.com/ecnerwala/cp-book/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-github.io-blue?logo=github)](https://ecnerwala.github.io/cp-book/)
[![coverage](https://img.shields.io/endpoint?url=https%3A%2F%2Fecnerwala.github.io%2Fcp-book%2Fcoverage%2Fbadge.json)](https://ecnerwala.github.io/cp-book/coverage/)

This is my library of reference code for competitive programming. The goal is to
write generic, fast, and clean algorithm implementations for use in contests
like CodeForces or ICPC.

The library lives in `src/` as standalone headers organized by area
(`#include "fft/series.hpp"`, `#include "ds/seg_tree.hpp"`, compiled with
`-I src`), all inside `namespace wala`. The library requires C++23 (the whole toolchain
builds with `-std=c++23`, e.g. g++ >= 13). Browsable source, verification
results, and coverage are hosted at https://ecnerwala.github.io/cp-book/.

## Building and testing

```sh
cmake -B build && cmake --build build
ctest --test-dir build          # or ./build/tests
```

## Bundling

Unfortunately, some files (e.g. full power series with all fft engines) are too long to fit in CF's submission limit,
so we need to bundle/minify them for submission.

```sh
scripts/bundle.py verify/fft/convolution_mod.test.cpp > submission.cpp
scripts/bundle.py fft/series.hpp ds/seg_tree.hpp   # bundle headers together
scripts/bundle.py --minify fft/series.hpp          # compiler-directed minification
scripts/bundle.py --all                            # pregenerate all headers
```

To run the bundler offline, first initialize the UV cache by running
```sh
uv sync --script scripts/bundle.py
```

Then, you can run
```sh
UV_OFFLINE=1 scripts/bundle.py ...
uv run --offline --script scripts/bundle.py ...
```

Upgrade the bundler's lockfile with
```sh
uv lock --script scripts/bundle.py --upgrade
```
(Make sure to rerun the sync afterwards.)

## Library Checker verification

`verify/` holds [Library Checker](https://judge.yosupo.jp/) solutions,
verified in CI with
[competitive-verifier](https://github.com/competitive-verifier/competitive-verifier).

```sh
# Run verification locally
uvx competitive-verifier oj-resolve --include src verify --exclude third_party \
    --config .competitive-verifier/config.toml > verify_files.json
uvx competitive-verifier verify --verify-json verify_files.json
```

## Contest tooling

`contest/` holds standalone tooling for contests: a problem-directory
template (Makefile with sanitizers and precompiled headers), a layered
`.template` instantiator (`make_prob.py`), and a
[Competitive Companion](https://github.com/jmerle/competitive-companion)
listener (`download_prob.py`). See [contest/README.md](contest/README.md).
Contest programs are standalone; vendor book code in as needed (e.g. with
`scripts/bundle.py`).

## License and attribution

All code is written by me and CC0 licensed unless otherwise noted in the
file. Inspiration is largely drawn from
[KACTL](https://github.com/kth-competitive-programming/kactl/) and other
references.
