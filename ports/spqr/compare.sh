#!/bin/bash
# usage: compare.sh <bin> <start_seed> <end_seed>
bin=$1; fails=0
for seed in $(seq $2 $3); do
  python3 gen.py $seed > /tmp/in_$$.txt
  ./cpp/dump_cpp < /tmp/in_$$.txt > /tmp/a_$$.txt
  $bin < /tmp/in_$$.txt > /tmp/b_$$.txt 2>/tmp/err_$$.txt
  if ! cmp -s /tmp/a_$$.txt /tmp/b_$$.txt; then echo "MISMATCH seed=$seed"; head -c 300 /tmp/err_$$.txt; fails=$((fails+1)); [ $fails -ge 3 ] && break; fi
done
echo "done $2..$3 fails=$fails"
