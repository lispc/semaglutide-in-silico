#!/bin/bash
# exp-C chain-length series queue: c12/c14/c16/c20/c22 diacid x 3 replicas x 100 ns.
#   GPU 0: c12 r1, c12 r2, c12 r3, c14 r1, c14 r2
#   GPU 1: c14 r3, c16 r1, c16 r2, c16 r3, c20 r1
#   GPU 2: c20 r2, c20 r3, c22 r1, c22 r2, c22 r3
# GPU 3 reserved for another task. 100 ns each (50,000,000 steps @ 2 fs).
# All GPUs verified free at launch (2026-07-20); no wait condition needed.
set -u
REPO=/home/scroll/personal/semaglutide-in-silico
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd "$REPO"
PY=/home/scroll/miniforge3/envs/gmx/bin/python3

echo "[$(date '+%F %T')] exp-C chain-series queue started (PID $$). GPU0/1/2, 5 jobs each."

run_rep() {  # system replica gpu
    local sys=$1 rep=$2 gpu=$3
    local outdir="exps/exp-C/md/${sys}/rep${rep}"
    mkdir -p "$outdir"
    echo "[$(date '+%F %T')] START ${sys} rep${rep} on GPU ${gpu}"
    CUDA_VISIBLE_DEVICES=$gpu "$PY" -u exps/exp-C/md/run_md.py \
        --system "$sys" --replica "$rep" --gpu 0 --nsteps 50000000 \
        > "$outdir/md_output.log" 2>&1
    echo "[$(date '+%F %T')] DONE  ${sys} rep${rep} on GPU ${gpu} (exit $?)"
    sleep 10
}

( run_rep c12_diacid 1 0; run_rep c12_diacid 2 0; run_rep c12_diacid 3 0; \
  run_rep c14_diacid 1 0; run_rep c14_diacid 2 0 ) &
( run_rep c14_diacid 3 1; run_rep c16_diacid 1 1; run_rep c16_diacid 2 1; \
  run_rep c16_diacid 3 1; run_rep c20_diacid 1 1 ) &
( run_rep c20_diacid 2 2; run_rep c20_diacid 3 2; run_rep c22_diacid 1 2; \
  run_rep c22_diacid 2 2; run_rep c22_diacid 3 2 ) &
wait
echo "[$(date '+%F %T')] exp-C chain-series queue complete: 15 jobs finished."
