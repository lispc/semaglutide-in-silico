#!/bin/bash
# Run inside screen/tmux for persistence
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-F/gmx
export CUDA_VISIBLE_DEVICES=2

# NPT equilibration (continue from checkpoint if exists)
if [ ! -f npt.gro ]; then
    echo "$(date): Starting NPT from checkpoint..."
    gmx mdrun -deffnm npt -ntmpi 1 -nb gpu -pme gpu -update gpu -cpi npt.cpt
    echo "$(date): NPT complete."
else
    echo "$(date): NPT already complete."
fi

# Production MD
if [ ! -f md.cpt ]; then
    echo "$(date): Starting production MD..."
    gmx grompp -f md.mdp -c npt.gro -p system.top -n system.ndx -o md.tpr -maxwarn 5
    gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -update gpu
else
    echo "$(date): Continuing production MD from checkpoint..."
    gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -update gpu -cpi md.cpt
fi

echo "$(date): All done!"
