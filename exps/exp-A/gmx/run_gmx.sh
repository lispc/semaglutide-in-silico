#!/bin/bash
# GROMACS MD for cross-validation. Usage: ./run_gmx.sh wt|aib8 GPU_ID
set -e
LABEL=$1
GPU=$2
DIR=$(dirname $0)/${LABEL}
cd $DIR

source ~/miniforge3/etc/profile.d/conda.sh
conda activate gmx

COORD="${LABEL}_coord.gro"

echo "=== ${LABEL} GPU ${GPU}: EM ==="
gmx grompp -f em.mdp -c ${COORD} -p ${LABEL}.top -o em.tpr -maxwarn 50
gmx mdrun -deffnm em -v -gpu_id ${GPU} -ntmpi 1 -ntomp 8 -pin off 2>&1 | tail -3

echo "=== ${LABEL} GPU ${GPU}: NVT (0→310K) ==="
gmx grompp -f nvt.mdp -c em.gro -p ${LABEL}.top -o nvt.tpr -maxwarn 50
gmx mdrun -deffnm nvt -v -gpu_id ${GPU} -ntmpi 1 -ntomp 8 -pin off 2>&1 | tail -3

echo "=== ${LABEL} GPU ${GPU}: NPT ==="
gmx grompp -f npt.mdp -c nvt.gro -p ${LABEL}.top -o npt.tpr -maxwarn 50
gmx mdrun -deffnm npt -v -gpu_id ${GPU} -ntmpi 1 -ntomp 8 -pin off 2>&1 | tail -3

echo "=== ${LABEL} GPU ${GPU}: Production MD ==="
gmx grompp -f md.mdp -c npt.gro -p ${LABEL}.top -o md.tpr -maxwarn 50
gmx mdrun -deffnm md -v -gpu_id ${GPU} -ntmpi 1 -ntomp 8 -pin off 2>&1 | tail -3

echo "=== ${LABEL} DONE ==="
