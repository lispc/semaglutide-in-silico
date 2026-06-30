#!/usr/bin/env python3
"""
Convert membrane system from Amber (OpenMM) to GROMACS for cross-validation.
Uses equilibrated coordinates from OpenMM equilibration.

Follows best-practice.md:
  - nstlist=40, rlist=1.1
  - lincs_iter=2, lincs_order=6
  - emtol <= 100
  - semi-isotropic pressure coupling for membrane
  - no CPU pinning for multi-GPU
"""
import os
import sys
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
GMX_DIR = f"{EXP_F}/gmx"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
# Use equilibrated PDB from OpenMM equilibration as starting structure
EQUIL_PDB = f"{EXP_F}/md/membrane_equil/equilibrated.pdb"

os.makedirs(GMX_DIR, exist_ok=True)

print("=" * 60)
print("Membrane System GROMACS Setup")
print(f"Prmtop: {PRMTOP}")
print(f"Coords: {EQUIL_PDB}")
print(f"Output: {GMX_DIR}")
print("=" * 60)

# --- 1. Load Amber topology with equilibrated coordinates ---
print("\n[1/4] Loading Amber topology + equilibrated coordinates...")
amber = pmd.load_file(PRMTOP, EQUIL_PDB)
print(f"  Atoms: {len(amber.atoms)}")
print(f"  Residues: {len(amber.residues)}")

# --- 2. Save as GROMACS top + gro ---
print("\n[2/4] Converting to GROMACS format...")
top_path = f"{GMX_DIR}/system.top"
gro_path = f"{GMX_DIR}/system.gro"
amber.save(top_path, format="gromacs", overwrite=True)
amber.save(gro_path, format="gromacs", overwrite=True)
print(f"  Top: {top_path}")
print(f"  Gro: {gro_path}")

# --- 3. Write MDP files ---
print("\n[3/4] Writing MDP files...")

# Energy minimization
mdp_em = """
; Energy minimization for membrane system
integrator      = steep
nsteps           = 5000
emtol           = 100.0
emstep          = 0.01
nstxout         = 0
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
lincs_iter      = 2
lincs_order     = 6
"""

# NVT heating: 100K -> 310K with annealing
mdp_nvt = """
; NVT heating with annealing
integrator      = md
dt              = 0.002
nsteps           = 50000       ; 100 ps
nstxout-compressed = 5000
nstlog          = 5000
nstenergy       = 5000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
; Temperature coupling
tcoupl          = v-rescale
tc-grps         = Protein Lipid Water_and_Ions
tau-t           = 0.1 0.1 0.1
ref-t           = 310 310 310
; Annealing: 100K -> 310K over 100ps
annealing       = single single single
annealing-npoints = 2 2 2
annealing-time  = 0 100 0 100 0 100
annealing-temp  = 100 310 100 310 100 310
; No pressure coupling
pcoupl          = no
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
lincs_iter      = 2
lincs_order     = 6
; Neighbor list
nstlist         = 40
rlist           = 1.1
"""

# NPT equilibration: semi-isotropic for membrane
mdp_npt = """
; NPT equilibration for membrane (semi-isotropic)
integrator      = md
dt              = 0.002
nsteps           = 500000      ; 1 ns
nstxout-compressed = 5000
nstlog          = 5000
nstenergy       = 5000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
; Temperature coupling
tcoupl          = v-rescale
tc-grps         = Protein Lipid Water_and_Ions
tau-t           = 0.1 0.1 0.1
ref-t           = 310 310 310
; Semi-isotropic pressure coupling (membrane standard)
pcoupl          = C-rescale
pcoupltype      = semiisotropic
tau-p           = 2.0
ref-p           = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
refcoord-scaling = com
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
lincs_iter      = 2
lincs_order     = 6
; Neighbor list
nstlist         = 40
rlist           = 1.1
"""

# Production MD: NPT semi-isotropic, 200 ns
mdp_md = """
; Production MD for membrane system (200 ns)
integrator      = md
dt              = 0.002
nsteps           = 100000000   ; 200 ns
nstxout-compressed = 5000      ; 10 ps (match OpenMM DCD interval)
nstlog          = 5000
nstenergy       = 5000
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0
; Temperature coupling
tcoupl          = v-rescale
tc-grps         = Protein Lipid Water_and_Ions
tau-t           = 0.1 0.1 0.1
ref-t           = 310 310 310
; Semi-isotropic pressure coupling
pcoupl          = C-rescale
pcoupltype      = semiisotropic
tau-p           = 2.0
ref-p           = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
refcoord-scaling = com
pbc             = xyz
constraints     = h-bonds
constraint-algorithm = LINCS
lincs_iter      = 2
lincs_order     = 6
; Neighbor list
nstlist         = 40
rlist           = 1.1
; Remove center-of-mass motion
comm-mode       = linear
comm-grps       = Protein Lipid Water_and_Ions
; Initial velocities (different seed from OpenMM)
gen_vel         = yes
gen_temp        = 310
gen_seed        = 42
"""

for name, content in [("em", mdp_em), ("nvt", mdp_nvt), ("npt", mdp_npt), ("md", mdp_md)]:
    path = f"{GMX_DIR}/{name}.mdp"
    with open(path, "w") as f:
        f.write(content.strip() + "\n")
    print(f"  {path}")

# --- 4. Write index generation script ---
print("\n[4/4] Writing index generation script...")
index_script = f"""#!/bin/bash
# Generate index file for temperature coupling groups
# Must be run after grompp with make_ndx

cd {GMX_DIR}

echo "Generating index file for tc-grps..."

# Create index with lipid group
gmx make_ndx -f system.gro -o system.ndx << 'EOF'
r PA PC OL CHL
name 19 Lipid
r WAT K+ CL-
name 20 Water_and_Ions
q
EOF

echo "Index file: system.ndx"
echo "Groups: Protein, Lipid, Water_and_Ions"
"""

with open(f"{GMX_DIR}/make_index.sh", "w") as f:
    f.write(index_script)
os.chmod(f"{GMX_DIR}/make_index.sh", 0o755)

# --- 5. Write run script ---
run_script = f"""#!/bin/bash
# GROMACS membrane system MD (cross-validation vs OpenMM)
# GPU: 2 (OpenMM uses GPU 1)
# Follows best-practice: no CPU pinning, nstlist=40, semi-isotropic pressure

cd {GMX_DIR}
export CUDA_VISIBLE_DEVICES=2

echo "=== Step 1: Generate index ==="
bash make_index.sh

echo "=== Step 2: Energy minimization ==="
gmx grompp -f em.mdp -c system.gro -p system.top -n system.ndx -o em.tpr -maxwarn 2
gmx mdrun -deffnm em -ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu

echo "=== Step 3: NVT heating (100K -> 310K) ==="
gmx grompp -f nvt.mdp -c em.gro -p system.top -n system.ndx -o nvt.tpr -maxwarn 2
gmx mdrun -deffnm nvt -ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu

echo "=== Step 4: NPT equilibration (1 ns) ==="
gmx grompp -f npt.mdp -c nvt.gro -p system.top -n system.ndx -o npt.tpr -maxwarn 2
gmx mdrun -deffnm npt -ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu

echo "=== Step 5: Production MD (200 ns) ==="
gmx grompp -f md.mdp -c npt.gro -p system.top -n system.ndx -o md.tpr -maxwarn 2
gmx mdrun -deffnm md -ntmpi 1 -nb gpu -pme gpu -bonded gpu -update gpu

echo "Done!"
"""

with open(f"{GMX_DIR}/run_all.sh", "w") as f:
    f.write(run_script)
os.chmod(f"{GMX_DIR}/run_all.sh", 0o755)

print(f"\nSetup complete. Files in {GMX_DIR}:")
for f in sorted(os.listdir(GMX_DIR)):
    print(f"  {f}")
print("\nNext: cd exps/exp-F/gmx && bash run_all.sh")
