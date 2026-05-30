#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd /home/scroll/personal/semaglutide-in-silico
rep=$1
gpu=$((rep - 1))
mkdir -p exps/exp-C/md/full_sema/rep${rep}
CUDA_VISIBLE_DEVICES=$gpu python -u exps/exp-C/md/run_sema.py --replica $rep --gpu 0 > exps/exp-C/md/full_sema/rep${rep}/md_output.log 2>&1 &
echo "Full sema rep $rep on GPU $gpu (PID $!)"
