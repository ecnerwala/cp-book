#!/bin/bash
# usage: compare_lean.sh <start_seed> <end_seed>
fails=0
for seed in $(seq $1 $2); do
  python3 gen.py $seed > /tmp/in_$$.txt
  ./cpp/dump_np_cpp < /tmp/in_$$.txt > /tmp/a_$$.txt
  ./lean/.lake/build/bin/dump_lean < /tmp/in_$$.txt > /tmp/b_$$.txt 2>/tmp/err_$$.txt
  if ! cmp -s /tmp/a_$$.txt /tmp/b_$$.txt; then echo "MISMATCH seed=$seed"; head -c 300 /tmp/err_$$.txt; fails=$((fails+1)); [ $fails -ge 3 ] && break; fi
done
echo "done $1..$2 fails=$fails"
