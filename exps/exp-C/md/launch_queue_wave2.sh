#!/bin/bash
# exp-C wave-2 queue: c16_monoacid (true C16 monoacid, liraglutide chain) x3,
# ONE replica appended per GPU 0/1/2, each starting ONLY after that GPU's
# wave-1 worker (launch_queue_chains.sh subshell) has exited.
# Does NOT touch the running launch_queue_chains.sh (editing a running bash
# script is unsafe -- bash reads it incrementally).
#
# Worker->GPU mapping captured 2026-07-20 (from /proc/<python>/environ):
#   1446459 -> GPU0 (c12 r1-3, c14 r1-2)
#   1446460 -> GPU1 (c14 r3, c16 r1-3, c20 r1)
#   1446461 -> GPU2 (c20 r2-3, c22 r1-3)
# If a worker PID is dead or its cmdline no longer matches launch_queue_chains
# (e.g. queue restarted), that replica starts immediately.
set -u
REPO=/home/scroll/personal/semaglutide-in-silico
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd "$REPO"
PY=/home/scroll/miniforge3/envs/gmx/bin/python3

W0=1446459; W1=1446460; W2=1446461

wait_worker() {  # $1 = wave-1 worker PID
    while [ -d "/proc/$1" ] \
          && tr '\0' ' ' < "/proc/$1/cmdline" 2>/dev/null | grep -q "launch_queue_chains"; do
        sleep 600
    done
}

echo "[$(date '+%F %T')] wave-2 queue started (PID $$). Waiting on wave-1 workers: $W0(GPU0) $W1(GPU1) $W2(GPU2)"

run_rep() {  # system replica gpu
    local sys=$1 rep=$2 gpu=$3
    local outdir="exps/exp-C/md/${sys}/rep${rep}"
    mkdir -p "$outdir"
    echo "[$(date '+%F %T')] START ${sys} rep${rep} on GPU ${gpu}"
    CUDA_VISIBLE_DEVICES=$gpu "$PY" -u exps/exp-C/md/run_md.py \
        --system "$sys" --replica "$rep" --gpu 0 --nsteps 50000000 \
        > "$outdir/md_output.log" 2>&1
    echo "[$(date '+%F %T')] DONE  ${sys} rep${rep} on GPU ${gpu} (exit $?)"
}

( wait_worker $W0; echo "[$(date '+%F %T')] GPU0 wave-1 done"; run_rep c16_monoacid 1 0 ) &
( wait_worker $W1; echo "[$(date '+%F %T')] GPU1 wave-1 done"; run_rep c16_monoacid 2 1 ) &
( wait_worker $W2; echo "[$(date '+%F %T')] GPU2 wave-1 done"; run_rep c16_monoacid 3 2 ) &
wait
echo "[$(date '+%F %T')] wave-2 queue complete: 3 jobs finished."
