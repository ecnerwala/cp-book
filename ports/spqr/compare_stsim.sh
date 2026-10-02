#!/bin/bash
# usage: compare_stsim.sh <start_seed> <end_seed>
# Checks the simulation relation `StSim` (Spqr/StSim.lean) at every finishEdge boundary of the walk;
# needs `lake build check_stsim`.
fails=0
for seed in $(seq $1 $2); do
  python3 gen.py $seed > /tmp/ss_$$.txt
  r=$(./lean/.lake/build/bin/check_stsim < /tmp/ss_$$.txt 2>&1)
  case "$r" in
    "checks "*" bad 0") ;;
    *) echo "MISMATCH seed=$seed: $r" | head -c 600; echo; fails=$((fails+1)); [ $fails -ge 3 ] && break;;
  esac
done
echo "done $1..$2 fails=$fails"
