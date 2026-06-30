#!/bin/bash
# GROMACS membrane system MD (cross-validation vs OpenMM)
# GPU: 2 (OpenMM membrane uses GPU 1)
# Follows best-practice: no CPU pinning, nstlist=40, semi-isotropic pressure

set -e
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-F/gmx
export CUDA_VISIBLE_DEVICES=2

echo "=== Step 1: Generate index ==="
if [ ! -f system.ndx ]; then
    bash make_index.sh
fi

echo "=== Step 2: Energy minimization ==="
if [ ! -f em.gro ]; then
    gmx grompp -f em.mdp -c system.gro -p system.top -n system.ndx -o em.tpr -maxwarn 5
    gmx mdrun -deffnm em -ntmpi 1
    echo "EM done. PE=$(grep 'Potential Energy' em.log | tail -1)"
fi

echo "=== Step 3: NVT heating (100K -> 310K, 100 ps) ==="
if [ ! -f nvt.gro ]; then
    gmx grompp -f nvt.mdp -c em.gro -p system.top -n system.ndx -o nvt.tpr -maxwarn 5
    gmx mdrun -deffnm nvt -ntmpi 1 -nb gpu -pme gpu -update gpu
    echo "NVT heating done."
fi

echo "=== Step 4: NPT equilibration (1 ns, semi-isotropic) ==="
if [ ! -f npt.gro ]; then
    gmx grompp -f npt.mdp -c nvt.gro -p system.top -n system.ndx -o npt.tpr -maxwarn 5
    gmx mdrun -deffnm npt -ntmpi 1 -nb gpu -pme gpu -update gpu
    echo "NPT equil done."
fi

echo "=== Step 5: Production MD (200 ns) ==="
if [ ! -f md.cpt ]; then
    # Start from scratch
    gmx grompp -f md.mdp -c npt.gro -p system.top -n system.ndx -o md.tpr -maxwarn 5
    gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -update gpu
else
    # Restart from checkpoint
    gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -update gpu -cpi md.cpt
fi

echo "Done!"
