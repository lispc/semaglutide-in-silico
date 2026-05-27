#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md
tleap -f build_c18v2.in
