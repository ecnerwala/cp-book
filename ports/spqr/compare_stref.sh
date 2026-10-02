#!/bin/bash
# usage: compare_stref.sh <start_seed> <end_seed>
# Checks that the walk's S / P / R child lists equal the restriction of the st-order reference
# (Spqr/StRef.lean); needs `lake build check_stref`.
fails=0
for seed in $(seq $1 $2); do
  python3 gen.py $seed > /tmp/sr_$$.txt
  r=$(./lean/.lake/build/bin/check_stref < /tmp/sr_$$.txt 2>&1)
  case "$r" in
    *"bad 0") ;;
    *) echo "MISMATCH seed=$seed: $r" | head -c 400; echo; fails=$((fails+1)); [ $fails -ge 3 ] && break;;
  esac
done
echo "done $1..$2 fails=$fails"
