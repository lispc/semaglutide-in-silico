#!/bin/bash
# Chain-series batch watcher: wait for rep1 batch exit 0 -> launch rep2 -> rep3.
set -u
BASE=/home/scroll/zzhang/semaglutide-in-silico/exps/exp-C/analysis/mmgbsa/chain-series
SYSTEMS="c12_diacid c14_diacid c16_diacid c20_diacid c22_diacid"

wait_batch() {
  local R=$1
  while true; do
    local done=0
    for s in $SYSTEMS; do
      if grep -q "exit: 0" $BASE/$s/rep$R/FINAL_run.log 2>/dev/null; then done=$((done+1)); fi
    done
    [ "$done" -eq 5 ] && break
    # bail if no MMPBSA processes remain but batch incomplete
    if [ "$(ps -eo args | grep -c '[M]MPBSA.py.MPI')" -eq 0 ]; then
      echo "[batchwatcher] ERROR: processes gone, done=$done/5 for rep$R" ; exit 1
    fi
    sleep 120
  done
}

echo "[batchwatcher] $(date) waiting rep1 batch"
wait_batch 1
echo "[batchwatcher] $(date) rep1 batch done; launching rep2"
cd $BASE
for s in $SYSTEMS; do setsid nohup ./run_rep.sh $s 2 14 > /dev/null 2>&1 < /dev/null & sleep 3; done
wait_batch 2
echo "[batchwatcher] $(date) rep2 batch done; launching rep3"
for s in $SYSTEMS; do setsid nohup ./run_rep.sh $s 3 14 > /dev/null 2>&1 < /dev/null & sleep 3; done
wait_batch 3
echo "[batchwatcher] $(date) ALL DONE"
