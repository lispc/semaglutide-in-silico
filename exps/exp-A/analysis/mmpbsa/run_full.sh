#!/bin/bash
# Full MM-GBSA/PBSA + decomp run for one exp-A system.
# Usage: run_full.sh <wt|aib8> <nranks>
set -u
SYS=$1
NP=${2:-16}
BASE=/home/scroll/zzhang/semaglutide-in-silico/exps/exp-A/analysis/mmpbsa
export AMBERHOME=/home/scroll/miniforge3/envs/cgas-md
export PATH=$AMBERHOME/bin:$PATH
export PYTHONPATH=$BASE/vendor
cd $BASE/$SYS

if [ "$SYS" = "wt" ]; then TRAJ=wt_dry_last100ns.nc; else TRAJ=aib8_dry_last100ns.nc; fi

/usr/bin/mpirun -np $NP MMPBSA.py.MPI -O -i mmpbsa.in \
  -cp com_dry.prmtop -rp rec_dry.prmtop -lp lig_dry.prmtop \
  -y $TRAJ \
  -o FINAL_MMPBSA.dat -do FINAL_DECOMP.dat \
  -eo FINAL_MMPBSA.csv -deo FINAL_DECOMP.csv \
  > mmpbsa_run.log 2>&1
echo "exit: $?" >> mmpbsa_run.log
