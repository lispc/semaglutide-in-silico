#!/bin/bash
# Usage: run_rep.sh <variant> <rep> <nranks> [input] [outprefix]
set -u
V=$1; R=$2; NP=${3:-10}; IN=${4:-mmpbsa.in}; OUT=${5:-FINAL}
BASE=/home/scroll/zzhang/semaglutide-in-silico/exps/exp-D/analysis/mmgbsa
export AMBERHOME=/home/scroll/miniforge3/envs/cgas-md
export PATH=$AMBERHOME/bin:$PATH
export PYTHONPATH=/home/scroll/zzhang/semaglutide-in-silico/exps/exp-A/analysis/mmpbsa/vendor
cd $BASE/$V/rep$R
/usr/bin/mpirun -np $NP MMPBSA.py.MPI -O -i ../$IN \
  -cp ../com_dry.prmtop -rp ../rec_dry.prmtop -lp ../lig_dry.prmtop \
  -y ${V}_rep${R}_dry.nc \
  -o ${OUT}_MMPBSA.dat -do ${OUT}_DECOMP.dat \
  -eo ${OUT}_MMPBSA.csv \
  > ${OUT}_run.log 2>&1
echo "exit: $?" >> ${OUT}_run.log
