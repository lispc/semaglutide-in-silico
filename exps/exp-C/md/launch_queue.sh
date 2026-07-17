#!/bin/bash
# exp-C production queue: HSA + free C18 monoacid/diacid at FA3.
#   GPU 2: c18_monoacid x 3 replicas (serial)
#   GPU 3: c18_diacid  x 3 replicas (serial)
#   100 ns each (50,000,000 steps @ 2 fs). GPU 0/1 reserved for exp-D, never used.
# The queue WAITS for exp-A MD (common/scripts/run_md.py) to finish first.
set -u
REPO=/home/scroll/personal/semaglutide-in-silico
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd "$REPO"
PY=/home/scroll/miniforge3/envs/gmx/bin/python3

echo "[$(date '+%F %T')] exp-C queue started (PID $$). Waiting for exp-A (common/scripts/run_md.py) to finish..."
while pgrep -f "common/scripts/run_md.py" >/dev/null; do
    sleep 600
done
echo "[$(date '+%F %T')] exp-A done. Starting exp-C jobs: GPU2=c18_monoacid x3, GPU3=c18_diacid x3."

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

( for rep in 1 2 3; do run_rep c18_monoacid "$rep" 2; done ) &
( for rep in 1 2 3; do run_rep c18_diacid  "$rep" 3; done ) &
wait
echo "[$(date '+%F %T')] exp-C queue complete: 6 jobs finished."
