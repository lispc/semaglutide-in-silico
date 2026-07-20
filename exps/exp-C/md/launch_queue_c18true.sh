#!/bin/bash
# exp-C c18true_diacid (true C18 diacid = HOOC-(CH2)16-COOH, 54 atoms) queue:
# 3 replicas x 100 ns on GPU 3. Waits for the exp-G pilot (any process
# matching "exp-G") to finish before touching GPU 3.
# This script's own cmdline contains "exp-C", never "exp-G" -> no self-match.
set -u
REPO=/home/scroll/personal/semaglutide-in-silico
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd "$REPO"
PY=/home/scroll/miniforge3/envs/gmx/bin/python3

echo "[$(date '+%F %T')] c18true queue started (PID $$). Waiting for exp-G pilot to finish..."
while pgrep -f "exp-G" >/dev/null; do
    sleep 600
done
echo "[$(date '+%F %T')] no exp-G process: GPU 3 free. Starting c18true_diacid x3."

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

for rep in 1 2 3; do
    run_rep c18true_diacid "$rep" 3
done
echo "[$(date '+%F %T')] c18true queue complete: 3 jobs finished."
