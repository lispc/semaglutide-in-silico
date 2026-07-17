#!/bin/bash
# exp-D production queue (restart 2026-07-17): 12 jobs = 4 variants x 3
# replicas, 100 ns each (50M steps). Serial per GPU:
#   GPU 0: no_linker rep1-3, gglu_1oeg rep1-3
#   GPU 1: gglu_2oeg rep1-3, gglu_3oeg rep1-3
# GPU 2/3 are busy with exp-A -- do not touch.
# Launch: setsid nohup bash launch_queue.sh >/dev/null 2>&1 &
PY=/home/scroll/miniforge3/envs/gmx/bin/python3
MD=/home/scroll/personal/semaglutide-in-silico/exps/exp-D/md

run_queue() {
  local gpu=$1; shift
  for spec in "$@"; do
    sys=${spec%%:*}; rep=${spec##*:}
    echo "[$(date '+%F %T')] GPU$gpu START $sys rep$rep"
    "$PY" "$MD/run_md.py" --system "$sys" --replica "$rep" --gpu "$gpu" --nsteps 50000000
    echo "[$(date '+%F %T')] GPU$gpu EXIT $? $sys rep$rep"
  done
  echo "[$(date '+%F %T')] GPU$gpu QUEUE DONE"
}

run_queue 0 no_linker:1 no_linker:2 no_linker:3 gglu_1oeg:1 gglu_1oeg:2 gglu_1oeg:3 >> "$MD/queue_gpu0.log" 2>&1 &
run_queue 1 gglu_2oeg:1 gglu_2oeg:2 gglu_2oeg:3 gglu_3oeg:1 gglu_3oeg:2 gglu_3oeg:3 >> "$MD/queue_gpu1.log" 2>&1 &
wait
