#!/bin/bash
# Auto-continue script: NVT -> NPT -> Production
# Runs in background, checks completion and launches next stage

source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate gmx
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-F/gmx
export CUDA_VISIBLE_DEVICES=2

# Wait for NVT to complete
while true; do
    if grep -q "Writing final coordinates" nvt.log 2>/dev/null; then
        echo "$(date): NVT complete. Starting NPT..."
        break
    fi
    if ! pgrep -f "mdrun.*nvt" > /dev/null; then
        echo "$(date): NVT process not found, checking if complete..."
        if [ -f nvt.gro ]; then
            echo "$(date): NVT gro exists, assuming complete."
            break
        fi
        echo "$(date): NVT failed or incomplete. Exiting."
        exit 1
    fi
    sleep 60
done

# NPT equilibration
if [ ! -f npt.gro ]; then
    gmx grompp -f npt.mdp -c nvt.gro -p system.top -n system.ndx -o npt.tpr -maxwarn 5
    gmx mdrun -deffnm npt -ntmpi 1 -nb gpu -pme gpu -update gpu
    echo "$(date): NPT complete."
fi

# Wait for NPT to complete
while true; do
    if grep -q "Writing final coordinates" npt.log 2>/dev/null; then
        echo "$(date): NPT complete. Starting production..."
        break
    fi
    if ! pgrep -f "mdrun.*npt" > /dev/null; then
        if [ -f npt.gro ]; then
            break
        fi
        echo "$(date): NPT failed. Exiting."
        exit 1
    fi
    sleep 60
done

# Production MD
if [ ! -f md.cpt ]; then
    gmx grompp -f md.mdp -c npt.gro -p system.top -n system.ndx -o md.tpr -maxwarn 5
    gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -update gpu
else
    gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -update gpu -cpi md.cpt
fi

echo "$(date): Production MD complete or checkpointed."
