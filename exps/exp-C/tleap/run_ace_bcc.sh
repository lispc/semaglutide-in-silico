#!/bin/bash
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate cgas-md
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap
rm -f ANTECHAMBER* ATOMTYPE* sqm.* 2>/dev/null
antechamber -i linker_ace.mol2 -fi mol2 -o linker_ace_bcc.mol2 -fo mol2 \
  -c bcc -at gaff2 -nc -1 -rn LFA -pf yes 2>&1 | tail -5
echo "Exit: $?"
ls -lh linker_ace_bcc.mol2 2>/dev/null
