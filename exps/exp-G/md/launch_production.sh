#!/bin/bash
# exp-G production queue: 3 replicas x 100 ns (50M steps each), GPU 3.
# rep1 starts immediately (shares GPU 3 with c18true queue until it drains).
# rep2/rep3 follow sequentially. Launch: setsid nohup bash launch_production.sh &
PY=/home/scroll/miniforge3/envs/gmx/bin/python3
MD=/home/scroll/personal/semaglutide-in-silico/exps/exp-G/md

for rep in 1 2 3; do
  seed=$((rep * 101))
  echo "[$(date '+%F %T')] START rep$rep (seed $seed)"
  CUDA_VISIBLE_DEVICES=3 "$PY" "$MD/run_pilot.py" --gpu 0 --nsteps 50000000 --out "rep$rep" --seed "$seed" >> "$MD/queue_production.log" 2>&1
  echo "[$(date '+%F %T')] EXIT $? rep$rep"
done
echo "[$(date '+%F %T')] PRODUCTION QUEUE DONE"
