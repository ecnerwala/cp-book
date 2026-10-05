#!/bin/bash
# usage: bench_all.sh [reps]   -- runs every bench binary on every bench_in/*.txt; prints "input | binary | line"
R=${1:-7}
cd "$(dirname "$0")"
for f in grid random sp tree random1m; do
  for b in cpp/bench_cpp rust/target/release/bench_rs zig/bench_zig zig/bench_zig_c zig/bench_zig_fast zig/bench_zig_fast_c; do
    ./$b $R < bench_in/$f.txt 2>/dev/null | sed "s@^@$f | $(basename $b) | @"
  done
done
