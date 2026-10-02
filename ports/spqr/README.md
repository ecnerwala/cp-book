# SPQR tree: Rust and Zig ports of `cp-book/src/graph/spqr_tree.hpp`

Five implementations of the same algorithm, all differential-tested to be byte-identical to the C++
on every output array (`vert_index`, `edge_index`, `par`, `subtree_end`, `types`, `orig_id`, `ch`,
`node_verts`, `vert_par_nv`, `node_edges`, `node_adj`, `node_planar`, `ne_rot_adj`, and
"non-planar build == planar build minus planarity"), on 300 random graphs per variant in debug (all asserts on)
plus 600 in release, covering disconnected graphs, isolated vertices, bridges,
self-loops, parallel edges, cycles, rigid graphs, `ternarize` on/off, and empty / full / prefix
`vert_order` / `edge_order`.

| file | style | ids | indexing / storage |
|---|---|---|---|
| `rust/src/spqr_tree.rs` | **faithful**: line-by-line port | `i32`, `-1` = none | `Vec`, checked indexing |
| `rust/src/spqr_tree_idiomatic.rs` | **idiomatic**: restructured as one would write it in Rust | `usize` locally, `Option<Idx>` (4-byte niche) stored | `Vec`, iterators, checked |
| `rust/src/spqr_tree_fast.rs` | **fast**: faithful port + unchecked indexing | `i32`, `-1` | `UVec`/`UCsr` (`get_unchecked`, asserted in debug) |
| `zig/spqr_tree.zig` | **faithful** | `i32`, `-1` | `ArrayList`, `Allocator` |
| `zig/spqr_tree_fast.zig` | **fast**: faithful + fixed-capacity stacks, vector fill, comptime layout | `i32`, `-1` | `Stack(T)` = `[]T` + `len`, preallocated |

## Layout

    rust/src/spqr_tree.rs             faithful Rust (SpqrTree / PlanarSpqrTree)
    rust/src/spqr_tree_idiomatic.rs   idiomatic Rust (same public build() shape, Option<Idx> / Csr<T> outputs)
    rust/src/spqr_tree_fast.rs        fast Rust (build / build_planar free functions, same output types as faithful)
    rust/src/dump.rs                  differential-test harness: dump_rs [faithful|fast|idiomatic] < graph
    rust/src/bench.rs                 benchmark harness (all three variants):  bench_rs [reps] < graph
    zig/spqr_tree.zig                 faithful Zig (Zig 0.16.0)
    zig/spqr_tree_fast.zig            fast Zig
    zig/dump.zig, zig/dump_fast.zig   differential-test harnesses
    zig/bench*.zig                    benchmark harnesses (bench_c* link libc and use std.heap.c_allocator)
    cpp/dump.cpp, cpp/bench.cpp       the same harnesses for the original C++ (reference output)
    gen.py, gen_big.py, gen_bench.py  random test / benchmark generators
    compare.sh <bin> <seed_from> <seed_to>   differential test driver;  bench_all.sh [reps]  runs everything

Build:

    cd rust && cargo build --release                     # target/release/{dump_rs,bench_rs}
    cd zig && zig build-exe dump.zig -O ReleaseFast && zig build-exe dump_fast.zig -O ReleaseFast
             zig build-exe bench_fast.zig -O ReleaseFast   # add -lc for bench_c*.zig
    g++ -std=c++23 -O2 -I ../../src cpp/dump.cpp -o cpp/dump_cpp

Differential test (each variant was run with seeds 0..300 in debug and 5000..5600 in release):

    ./compare.sh "./rust/target/release/dump_rs idiomatic" 5000 5600
    ./compare.sh ./zig/dump_zig_fast_release 5000 5600

## API

C++:
```cpp
auto t = wala::spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
auto p = wala::planar_spqr_tree::build(NV, edges, ternarize, vert_order, edge_order);
```
Rust (faithful and fast share the output types; fast exposes `spqr_tree_fast::{build, build_planar}`):
```rust
let t = SpqrTree::build(nv, &edges, ternarize, &vert_order, &edge_order);   // edges: &[[i32; 2]]
let p = PlanarSpqrTree::build(nv, &edges, ternarize, &vert_order, &edge_order);
p.node_planar, p.ne_rot_adj, and p.par / p.ch / ... via Deref<Target = SpqrTree>
```
Rust idiomatic (`usize` inputs, `Option<Idx>` where the C++ has `-1`, typed `NodeType` / `Csr<T>` outputs):
```rust
let t = spqr_tree_idiomatic::SpqrTree::build(nv, &edges, ternarize, &vert_order, &edge_order);  // edges: &[[usize; 2]]
t.par[i]: Option<Idx>,  t.types[i]: NodeType,  t.ch: Csr<Idx>,  t.node_verts: Csr<NodeVert { vert: Idx, par_nv: Option<Idx> }>
```
Zig (both variants):
```zig
var t = try SpqrTree.build(gpa, NV, edges, ternarize, vert_order, edge_order);   // edges: []const [2]i32
defer t.deinit(gpa);
var p = try PlanarSpqrTree.build(gpa, NV, edges, ternarize, vert_order, edge_order);
defer p.deinit(gpa);
```

## Faithful ports (where they deviate in *form*, never in output)

* Ids stay `i32` (with `-1` = none), and `csr<T>` stays `{bounds, dat}`.
* The C++ `[&]` lambdas capturing ~30 locals became methods on three state structs
  (`LowvalDfs`, `Builder`, `Relabel`), one per phase, with the same names.
* `template <bool with_planarity> build_impl` is a const generic `Builder<const WP: bool>`
  in Rust and `fn Builder(comptime WP: bool) type` in Zig.
* `std::expected<tstack_planarity, tstack_nonplanarity>` is `Result<_, TstackNonplanarity>` in Rust
  and `?TstackPlanarity` in Zig.
* Rust is allocation-free beyond `Vec`s; Zig takes an `Allocator`, returns `Allocator.Error!`, and the
  result owns its slices (`deinit`).
* Two Zig pitfalls that bit during porting (commented in the source): result-location aliasing
  (`lowvals.* = .{ a, f(lowvals[0]) }` is not a copy), and `&(opt orelse ..)` pointing at a temporary.

## Idiomatic Rust (`spqr_tree_idiomatic.rs`)

Same three phases and identical output, restructured rather than transliterated:

* **`Idx` = `-1` as a niche.** `struct Idx(NonZeroU32)` stores `!i`, so `u32::MAX` (the C++ `-1`) is
  the forbidden value and `Option<Idx>` is 4 bytes with `None` bit-identical to `-1`. `Ord` is
  implemented on the decoded value. Local arithmetic is in `usize`; stored arrays use `Idx` /
  `Option<Idx>`, so every `-1` sentinel and `x == -1` test became an `Option`.
* **State split instead of one 30-field struct.** The captured locals cluster into `Items`
  (item arrays + `alloc` / `concat`), `TStacks<WP>` (the triangle stacks + `merge_tops` / `flip`),
  `Planarity` (quarter-edge matches, per-tstack planarity shadow stack) and per-phase drivers
  (`LowvalDfs`, `Builder<WP>`, `Relabel<WP>`) that compose them with disjoint borrows.
* **Typed instead of packed.** `NodeType` enum with `is_node()`, `Key` (the `3 * (lowval + 2) + kind`
  edge key) with `encode` / `decode`, `Span` / `FlipItem` for the bit-packed tstack entries,
  `Result<TstackPlanarity, Nonplanar>` for `std::expected`.
* **Iterators for the linear passes.** `in_order(n, order)` replaces the `for_each_in_order` callback,
  `Csr::bucket(n, items, key)` is the counting sort used for `by_key` / `by_src` / adjacency,
  `prefix_sums` builds CSR bounds. The DFS / walk loops stay explicit loops over small frame structs.
* `expect("...")` documents each invariant the C++ relied on silently (e.g. "in a dfs", "child pending").
  No `unsafe`, clippy-clean.
* Frames and out-edges use `Idx` / `Range<u32>` fields so their sizes match the C++ structs (this was
  measurable: `usize` everywhere cost ~5–10%).

## Fast Rust (`spqr_tree_fast.rs`)

The faithful port with two mechanical changes; the algorithm text is otherwise identical:

* All hot-path arrays are `UVec<T>` / `UCsr<T>`: `Index`/`IndexMut` are `get_unchecked` in release
  and `debug_assert!`-checked in debug; `top()` / `pop_u()` likewise. This is the entire `unsafe`
  surface (5 sites, all in the `UVec` impl). Every index is one the C++ computes and uses unchecked,
  so the invariant is exactly the C++'s ("`build` never indexes out of range"), and the debug build
  is the checked version of it.
* Bounds checks the compiler *can* remove stay safe: `usize` loop indices, iterators over `edges`
  / `order`, single `assert!`s on input lengths. The DFS / tstack random accesses (`item_vs[ch_nxt[i]]`,
  `tstack[top_depth]`, ...) are the ones that need `UVec`.
* `Tstack` no longer carries the planarity payload; a parallel `tstack_planarity` vector holds it
  (empty when `WP == false`), so the non-planar build's tstack entries are 28 bytes instead of 80.

## Fast Zig (`spqr_tree_fast.zig`)

The faithful port with the three fixes the profile pointed at (see below):

* `Stack(T)` = `{ buf: []T, len: usize }` over a slice preallocated to the exact upper bound the C++
  `reserve`s (`NV`, `NE`, `1 + NV + 2 NE`, `NV + NE`, ...). Pushes are `assert(len < buf.len)` +
  store, so they inline and never allocate; `items()` / `pop()` / `top()` are the ArrayList
  operations the code used. The bound assertions run in debug; the capacities are the same formulas
  the C++ uses for `reserve`, which it relies on for pointer stability.
* `fill(slice, val)` writes 8-wide vectors (with an `asm volatile` barrier so LLVM's loop-idiom pass
  doesn't turn it back into `memset`) instead of `@memset`, whose `compiler_rt` implementation is
  ~4.6x slower than glibc's here for non-zero patterns.
* `planarity: if (WP) TstackMaybePlanarity else void` on the tstack entry: the non-planar build's
  entries shrink from 80 to 28 bytes (the Zig equivalent of the C++ `std::conditional_t<..., monostate>`).
* Allocator is still the caller's; `bench_c*.zig` show `std.heap.c_allocator` (glibc malloc), which is
  worth ~5–10% over the default `smp_allocator` for this allocation pattern.

## Performance

Release builds, best of 7, build time only (no I/O); Xeon 8559C; g++ 15.2 `-O2`; rustc 1.97 `-O`
(fat LTO / codegen-units=1 made no difference); zig 0.16 `ReleaseFast`. `+c` = `-lc` and
`std.heap.c_allocator` instead of the default `smp_allocator`. Times are `spqr / planar` ms.

    graph                    C++   Rust faithful      Rust idiom       Rust fast    Zig faithful  Zig faithful+c        Zig fast      Zig fast+c
    grid             62 /   89 ms     80 /  101 ms    105 /  134 ms     73 /   95 ms     97 /  128 ms     87 /  110 ms     77 /   95 ms     69 /   89 ms
    random          125 /  151 ms    136 /  162 ms    181 /  215 ms    129 /  151 ms    166 /  193 ms    149 /  181 ms    131 /  155 ms    126 /  153 ms
    series-par      126 /  161 ms    139 /  179 ms    176 /  223 ms    131 /  170 ms    182 /  223 ms    164 /  204 ms    130 /  171 ms    127 /  165 ms
    tree+10%         73 /   85 ms     80 /   93 ms    110 /  113 ms     76 /   86 ms     97 /  108 ms     87 /  101 ms     75 /   87 ms     71 /   88 ms
    random 1M       874 / 1024 ms    965 / 1142 ms   1258 / 1581 ms    918 / 1026 ms   1086 / 1241 ms    992 / 1176 ms    902 / 1043 ms    860 / 1015 ms

Summary: fast Rust and fast Zig (+c_allocator) are at parity with the C++ (within ~5%) on the
memory-bound graphs (random, series-parallel, tree, 1M); on the small-working-set grid the C++ keeps a
~15% lead over both. Faithful Rust is 5–20% behind, faithful Zig 25–35%, idiomatic Rust ~30% behind
faithful Rust.

GCC vs LLVM: the same C++ built with clang++-19 `-O2` gives grid 70 / 93, random 146 / 168, tree 80 / 94 ms
(g++: 62 / 79, 123 / 145, 74 / 86; fast Rust: 71 / 91, 132 / 156, 73 / 86), i.e. fast Rust matches or beats
LLVM-compiled C++ everywhere, and the residual gap to g++ on grid is GCC's backend, not the language.

Graphs: random multigraph (rigid), grid (planar rigid), series-parallel (S/P heavy), tree + 10% extra
edges, all NV≈200k / NE≈400k, and random NV=1M / NE=2M. Generated by `gen_bench.py`.

### Profiling notes (perf)

* Faithful Zig's ~30% gap vs C++ was: `compiler_rt.memset` (7–9%), `ArrayList.append` not inlining
  (~13%; `std::vector::push_back` is), and the default allocator (~10%). All three are addressed in
  `spqr_tree_fast.zig` (+ `c_allocator`).
* Faithful Rust's 5–20% gap is bounds checks (`build_impl<true>` has 235 panic call sites, `push_item`
  71 on 1375 instructions); `spqr_tree_fast.rs` removes them and is at parity with C++ on the
  memory-bound graphs.
* Idiomatic Rust is ~30% slower than the faithful port. Per-phase timing (grid): phase 1 (lowlink DFS +
  bucket sorts) 25 vs 18 ms, phase 2 (walk) 17/31 vs 13/22 ms, phase 3 (relabel) 56/69 vs 45/57 ms.
  The remaining cost is spread across `Option`/`expect` branches, `usize` CSR bounds (8 bytes vs 4
  on the heavily-touched `bounds` arrays), and `Idx` conversions; none of it is a single hot spot.
  Compacting the frame structs (`Idx` / `Range<u32>`) and removing a per-node `Vec` allocation in
  `push_item` recovered ~10%.
* The C++ `[&]` lambdas that g++ outlines (`push_item`, `finish_edge`, ...) reload captured scalars
  through the closure after every store (no `restrict`); Rust's `&mut self` with distinct fields
  doesn't have to. This is consistent with fast Rust matching C++ despite the `i32 -> usize` casts.
