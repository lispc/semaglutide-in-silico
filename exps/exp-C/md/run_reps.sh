#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd /home/scroll/personal/semaglutide-in-silico

for rep in 1 2 3; do
    gpu=$((rep - 1))
    echo "Launching rep $rep on GPU $gpu"
    CUDA_VISIBLE_DEVICES=$gpu nohup python exps/exp-C/md/run_md.py \
        --system c18_monoacid --replica $rep --gpu 0 \
        > exps/exp-C/md/c18_monoacid/rep${rep}/md_output.log 2>&1 &
done
echo "All 3 replicas launched"
