# AGENTS.md

Guidance for AI agents working in this repository (ecnerwala's competitive
programming reference library).

## Repository layout

- `src/`: the library — header-only C++ (`.hpp`), plus Catch2 unit tests
  (`src/*.test.cpp`). Headers are included with bare names (`#include
  "fft.hpp"`), resolved via `-I src`.
- `verify/`: Library Checker (https://judge.yosupo.jp/) solutions, one
  `<problem_slug>.test.cpp` per problem, verified by competitive-verifier.
- `scripts/bundle.py`: inlines library headers to produce a single
  submittable file (the docs site's bundled views come from oj-resolve's
  builtin bundling).
- `.competitive-verifier/config.toml`: compiler settings for verification
  (g++, `-std=c++26`, `-I src`, `read_macros = false`).

## Build, test, verify

```sh
# Unit tests (requires CMake >= 3.27 and a C++26 compiler, e.g. g++ >= 14)
cmake -B build && cmake --build build -j && ctest --test-dir build

# Library Checker verification (note: oj-resolve only sees git-tracked files,
# so `git add` new verify files first). The second (sanitizer) environment
# compiles through the scripts/toolchain/g++-sanitizer wrapper so it gets a
# distinct name in the results.
uvx competitive-verifier oj-resolve --include src verify --exclude third_party \
    --config .competitive-verifier/config.toml > verify_files.json
uvx competitive-verifier verify --verify-json verify_files.json

# Single-file submission with headers inlined
scripts/bundle.py verify/convolution_mod.test.cpp > submission.cpp
```

Compile a single file: `g++ -std=c++26 -O2 -I src file.cpp`.

## Style conventions

Match the existing code exactly. In particular:

- **Indent with tabs.**
- **No `using namespace std;`** — qualify with `std::` everywhere, including
  in `main`. Library names (e.g. `modnum`, `SuffixArray`) live in the global
  namespace or their own namespaces; that's intentional.
- **verify/ files** follow this template:

  ```cpp
  // competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/<slug>

  #include <bits/stdc++.h>

  #include "some_header.hpp"

  int main() {
  	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

  	// ...

  	return 0;
  }
  ```

  - Start with `#include <bits/stdc++.h>`, then a blank line, then the
    library includes (bare names, never `../src/...`).
  - Use the `// competitive-verifier: PROBLEM <url>` comment form, not
    `#define PROBLEM`.
- **Capitalize input variables** (`N`, `Q`, `A`, `B`, `S`, `F`, `G`...);
  loop indices and scratch stay lowercase.
- **Read input directly into library types** (`poly<E>`, `power_series`,
  `dirichlet_series_prefix::st[]`) rather than into an intermediate
  `std::vector` first. If a type is missing a small piece of interface that
  would make this natural (a sized constructor, mutable iteration), add it to
  the header rather than working around it.
- Compact declare-and-read style: `int N, Q; std::cin >> N >> Q;` on one
  line; `std::vector<int> A(N); for (auto& x : A) std::cin >> x;`.
- Output with `" \n"[i+1==N]` for space-separated lines.
- Dispatch on ops explicitly: `if (op == 0) { ... } else if (op == 1) { ... }
  else assert(false);`.
- Prefer library idioms over raw loops when the library provides them (e.g.
  `seg_tree::point a(N-1); a >= 1; a--` with `a.c(0)` / `a.c(1)` children).
- Keep `main` minimal — all real logic belongs in `src/` headers.

## Library notes

- `src/fft.hpp`: everything lives in `namespace ecnerwala` (engines in
  `ecnerwala::fft`). `poly<E>` stores coefficients reversed but iterates and
  indexes in natural (x^0-first) order; `poly_evaluate` / `poly_interpolate`
  are in `ecnerwala`, not `ecnerwala::fft`.
- `src/dirichlet_series.hpp`: `div_vector_layout` is passed as a template
  non-type reference parameter, so declare it `static` and assign the size at
  runtime.
- Unit tests (`src/*.test.cpp`) use Catch2 and are built by CMake; they are
  not Library Checker verifications (they have no PROBLEM attribute).

## When adding a verify solution

1. Create `verify/<problem_slug>.test.cpp` following the template above,
   using book headers (don't vendor code into the file). If a problem can't
   be solved with book code, flag it instead of adding standalone code.
2. `git add` it, then run the verification commands above; all testcases
   must be AC locally before pushing.
3. Keep the CI workflow (`.github/workflows/ci.yml`) and README in sync
   if commands change.
