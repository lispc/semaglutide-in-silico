#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap
tleap -f build_mono.in
