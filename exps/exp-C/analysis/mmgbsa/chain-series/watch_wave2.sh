#!/bin/bash
set -u
BASE=/home/scroll/zzhang/semaglutide-in-silico/exps/exp-C/analysis/mmgbsa/chain-series
echo "[wave2] waiting 5 cases (c16m rep1-3, c18true rep1-2)"
while true; do
  n=0
  for f in $BASE/c16_monoacid/rep1/FINAL_run.log $BASE/c16_monoacid/rep2/FINAL_run.log $BASE/c16_monoacid/rep3/FINAL_run.log $BASE/c18true_diacid/rep1/FINAL_run.log $BASE/c18true_diacid/rep2/FINAL_run.log; do
    grep -q "exit: 0" $f 2>/dev/null && n=$((n+1))
  done
  [ "$n" -eq 5 ] && break
  if [ "$(ps -eo args | grep -c '[M]MPBSA.py.MPI')" -eq 0 ]; then echo "[wave2] processes gone, n=$n/5"; exit 1; fi
  sleep 180
done
echo "[wave2] $(date) ALL 5 DONE"
