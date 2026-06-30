#!/bin/bash
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-F/gmx
/home/scroll/miniforge3/envs/gmx/bin/python analysis/comprehensive_analysis.py > analysis/run.log 2>&1
echo "EXIT_CODE=$?" >> analysis/run.log
