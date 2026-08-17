# AGENTS.md

Guidance for AI agents working in this repository (ecnerwala's competitive
programming reference library).

## Repository layout

- `src/`: the library — header-only C++ (`.hpp`), plus Catch2 unit tests
  (`src/**/*.test.cpp`). Everything lives in `namespace wala`. Headers are
  organized into area folders (`num/`, `nt/`, `linalg/`, `ds/`, `seq/`,
  `tree/`, `graph/`, `geometry/`, `fft/`, `combo_games/`, plus
  uncategorized headers at the top level) and included with src-relative
  paths (`#include "fft/series.hpp"`, `#include "ds/seg_tree.hpp"`),
  resolved via `-I src`. Same-directory includes may use bare names.
- `verify/`: Library Checker (https://judge.yosupo.jp/) solutions, one
  `<problem_slug>.test.cpp` per problem, verified by competitive-verifier.
  Organized into the same area folders as `src/`
  (`verify/fft/convolution_mod.test.cpp`), chosen by the solution's primary
  header.
- `scripts/bundle.py`: inlines library headers to produce a single
  submittable file; takes one or more headers/solutions, `--minify` for a
  compiler-directed minification pass, and `--all` to pregenerate
  bundled+minified copies of every header (the docs site's bundled/minified
  views come from oj-resolve's builtin bundling).
- `.competitive-verifier/config.toml`: compiler settings for verification
  (g++, `-std=c++23`, `-I src`, `read_macros = false`).

## Build, test, verify

```sh
# Unit tests (requires CMake >= 3.27 and a C++23 compiler, e.g. g++ >= 13)
cmake -B build && cmake --build build -j && ctest --test-dir build

# Library Checker verification (note: oj-resolve only sees git-tracked files,
# so `git add` new verify files first). The second (sanitizer) environment
# compiles through the scripts/toolchain/g++-sanitizer wrapper so it gets a
# distinct name in the results.
uvx competitive-verifier oj-resolve --include src verify --exclude third_party \
    --config .competitive-verifier/config.toml > verify_files.json
uvx competitive-verifier verify --verify-json verify_files.json

# Single-file submission with headers inlined
scripts/bundle.py verify/fft/convolution_mod.test.cpp > submission.cpp
```

Compile a single file: `g++ -std=c++23 -O2 -I src file.cpp`.

## Style conventions

Match the existing code exactly. In particular:

- **Indent with tabs.**
- **No `using namespace`** — qualify with `std::` everywhere, including in
  `main`. The library lives in `namespace wala`; qualify library names with
  `wala::` in tests and verify solutions (e.g. `wala::modnum`,
  `wala::SuffixArray`, `wala::seg_tree::point`).
- **verify/ files** follow this template:

  ```cpp
  // competitive-verifier: PROBLEM https://judge.yosupo.jp/problem/<slug>

  #include <bits/stdc++.h>
  #include <cassert>

  #include "some_header.hpp"

  int main() {
  	std::ios_base::sync_with_stdio(false), std::cin.tie(nullptr);

  	// ...

  	return 0;
  }
  ```

  - Start with `#include <bits/stdc++.h>`, then `#include <cassert>`, then a
    blank line, then the library includes (src-relative paths like
    `"ds/seg_tree.hpp"`, never `../src/...`).
  - Always include `<cassert>` right after `<bits/stdc++.h>`, for
    consistency: it is not included in `<bits/stdc++.h>` in recent versions
    of g++, and ecnerwala loves `assert()`s.
  - Use the `// competitive-verifier: PROBLEM <url>` comment form, not
    `#define PROBLEM`.
- **Capitalize input variables** (`N`, `Q`, `A`, `B`, `S`, `F`, `G`...);
  loop indices and scratch stay lowercase.
- **Read input directly into library types** (`poly::vec<E>`, `series::vec`,
  `dirichlet_series::prefix::st[]`) rather than into an intermediate
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
- **Never skip braces** on loops/branches, with one exception: the body may
  be brace-less when it sits on the same line as the loop/branch and is the
  entire rest of that line (e.g.
  `for (int i = 0; i < n; i++) out[i] = a[i] + b[i];`).
- **Hanging indents use `(` then a newline**: when a declaration or call
  doesn't fit on one line, break immediately after the opening paren, put
  the arguments on following lines indented one extra tab, and close with
  `)` back at the outer indentation level. This keeps all indentation at
  whole multiples of tab (no aligning to the paren column):

  ```cpp
  template <fft::engine E>
  std::vector<typename E::value_type> multipoint(
  	const poly::vec<E>& p,
  	std::span<const typename E::value_type> pts
  ) {
  ```

## Comment style

Comments exist for the reader of the code as it stands — to understand it and
use it correctly. Concretely:

- **Describe the code as-is, never the change.** When editing, don't write
  comments that narrate what changed or contrast with previous behavior; that
  context belongs in the commit/PR message. (Occasionally a change embodies an
  important design decision worth recording — record the decision, not the
  diff.) Most implementation details are self-evident from the code and need
  no comment at all.
- **State a convention once, where it's established**, then rely on it —
  don't re-explain it at every use site. E.g. the bit-reversed transform
  convention and the span length rules are stated in the `fft_core` /
  `engine` preambles, and later functions just use them.
- **Contracts over prose.** Prefer short statements of what a function
  computes and requires (`Requires f.len() == N.`, allowed span lengths,
  aliasing rules) to paragraphs of explanation. Derivations/proofs, if worth
  keeping, go inside the function body next to the relevant code.
- **Start each sentence/clause on a new line** — this makes comments easier
  to read and edit in a line-based editor. Don't rewrap them into paragraphs.
- Plain inline math notation: `R[[x]]`, `2^k`, `omega`, `<P, S> = [x^0]
  P(1/x) S(x)`, half-open ranges `[Ofs[c], Ofs[c+1])`, python-style slices
  `prod[Ofs[c]:Ofs[c+1]]`. No LaTeX ceremony outside the file header.
- First-person-plural voice for design narration is fine ("We take the
  bit-reverse convention...", "We will store Q as a packed buffer...").
- Open design questions are recorded inline as `// TODO:` next to the code
  they concern. Leave existing TODOs alone unless you're resolving them.
- Avoid hand-aligned comment columns/tables unless the alignment clearly
  aids readability (e.g. the layer table in the fft.hpp header); aligned
  blocks are costly to maintain. No trailing whitespace.

## Library notes

- `src/fft/`: transform machinery in `wala::fft`, engines in
  `wala::fft::engines`, value types in `wala::series` / `wala::poly`.
  `poly::vec<E>` stores coefficients reversed but iterates and
  indexes in natural (x^0-first) order; `poly::multipoint` / `poly::interpolate`
  live in `wala::poly`.
- `src/nt/dirichlet_series.hpp`: `div_vector_layout` is passed as a template
  non-type reference parameter, so declare it `static` and assign the size at
  runtime.
- Unit tests (`src/*.test.cpp`) use Catch2 and are built by CMake; they are
  not Library Checker verifications (they have no PROBLEM attribute).

## When adding a verify solution

1. Create `verify/<area>/<problem_slug>.test.cpp` following the template above,
   using book headers (don't vendor code into the file). If a problem can't
   be solved with book code, flag it instead of adding standalone code.
2. `git add` it, then run the verification commands above; all testcases
   must be AC locally before pushing.
3. Keep the CI workflow (`.github/workflows/ci.yml`) and README in sync
   if commands change.
