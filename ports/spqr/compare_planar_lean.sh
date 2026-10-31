#!/bin/bash
# usage: compare_planar_lean.sh <start_seed> <end_seed>
# Compares dump_planar_lean against cpp/dump.cpp (+ dump_embed.cpp for the glued embedding) byte for byte.
fails=0
for seed in $(seq $1 $2); do
  python3 gen.py $seed > /tmp/in_$$.txt
  { ./cpp/dump_cpp < /tmp/in_$$.txt; ./cpp/dump_embed_cpp < /tmp/in_$$.txt; } > /tmp/a_$$.txt
  SPQR_EMBED=1 ./lean/.lake/build/bin/dump_planar_lean < /tmp/in_$$.txt > /tmp/b_$$.txt 2>/tmp/err_$$.txt
  if ! cmp -s /tmp/a_$$.txt /tmp/b_$$.txt; then echo "MISMATCH seed=$seed"; head -c 300 /tmp/err_$$.txt; fails=$((fails+1)); [ $fails -ge 3 ] && break; fi
done
echo "done $1..$2 fails=$fails"
# Lean-side specification check (Spqr.Planar) of the local and glued embeddings.
checkfails=0
for seed in $(seq $1 $2); do
  python3 gen.py $seed > /tmp/in_$$.txt
  if ! ./lean/.lake/build/bin/check_planar_lean < /tmp/in_$$.txt > /tmp/c_$$.txt 2>&1; then echo "CHECK FAIL seed=$seed"; head -c 300 /tmp/c_$$.txt; checkfails=$((checkfails+1)); [ $checkfails -ge 3 ] && break; fi
done
echo "check $1..$2 fails=$checkfails"
