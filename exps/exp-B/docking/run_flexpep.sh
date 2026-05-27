#!/bin/bash
# Rosetta FlexPepDocking refinement for exp-B
# Refines the positioned ECD-peptide complex
source /home/scroll/miniforge3/etc/profile.d/conda.sh
conda activate rosetta

INPUT="../structures/4ZGM_positioned.pdb"
OUTDIR="."

# Run FlexPepDocking in refinement mode
# Skip low-res preoptimization (best-practice 13: don't trust ab-initio for long peptides)
FlexPepDocking \
  -s $INPUT \
  -flexPepDocking:receptor_chain A \
  -flexPepDocking:peptide_chain B \
  -flexPepDocking:pep_refine \
  -flexPepDocking:lowres_preoptimize false \
  -nstruct 10 \
  -out:prefix flexpep_ \
  -out:path $OUTDIR \
  -mute all

echo "FlexPepDocking complete. Check $OUTDIR/flexpep_*.pdb"
