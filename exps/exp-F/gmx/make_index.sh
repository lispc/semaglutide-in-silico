#!/bin/bash
# Generate index file for temperature coupling groups
cd /home/scroll/personal/semaglutide-in-silico/exps/exp-F/gmx

echo "Generating index file for tc-grps..."

# Create index with lipid group
gmx make_ndx -f system.gro -o system.ndx << 'INDEXEOF'
r PA PC OL CHL
name 19 Lipid
r HOH K+ Cl-
name 20 Water_and_Ions
q
INDEXEOF

echo "Index file: system.ndx"
echo "Groups: Protein, Lipid, Water_and_Ions"
