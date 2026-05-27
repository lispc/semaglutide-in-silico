#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd /home/scroll/personal/semaglutide-in-silico
CUDA_VISIBLE_DEVICES=0 python exps/exp-C/md/run_md.py --system c18_monoacid --gpu 0 --nsteps 50000000
