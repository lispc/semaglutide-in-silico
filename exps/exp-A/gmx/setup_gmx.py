#!/usr/bin/env python3
"""Convert Amber prmtop/inpcrd to GROMACS and write mdp for cross-validation."""
import parmed as pmd, os

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_A = f"{REPO}/exps/exp-A"

for label in ["wt", "aib8"]:
    gmx_dir = f"{EXP_A}/gmx/{label}"
    os.makedirs(gmx_dir, exist_ok=True)

    # Load Amber topology
    if label == "wt":
        prmtop = f"{EXP_A}/tleap/wt.prmtop"
        inpcrd = f"{EXP_A}/tleap/wt.inpcrd"
    else:
        prmtop = f"{EXP_A}/md/aib8/aib8_modified.prmtop"
        inpcrd = f"{EXP_A}/md/aib8/aib8_modified.inpcrd"

    amber = pmd.load_file(prmtop, inpcrd)

    # Save as GROMACS
    amber.save(f"{gmx_dir}/{label}.top", format="gromacs", overwrite=True)
    amber.save(f"{gmx_dir}/{label}.gro", format="gromacs", overwrite=True)
    print(f"{label}: {len(amber.atoms)} atoms → {gmx_dir}/{label}.top, {label}.gro")

# Write GROMACS MDP files
mdp_em = """
integrator  = steep
nsteps      = 5000
emtol       = 1000
emstep      = 0.01
nstxout     = 0
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
pbc         = xyz
"""

mdp_nvt = """
integrator  = md
dt          = 0.002
nsteps      = 25000
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 5000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1 0.1
ref-t       = 100 100
pcoupl      = no
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
"""

mdp_npt = """
integrator  = md
dt          = 0.002
nsteps      = 150000
nstxout-compressed = 5000
nstlog      = 5000
nstenergy   = 5000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1 0.1
ref-t       = 310 310
pcoupl      = C-rescale
pcoupltype  = isotropic
tau-p       = 2.0
ref-p       = 1.0
compressibility = 4.5e-5
refcoord-scaling = com
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
"""

mdp_md = """
integrator  = md
dt          = 0.002
nsteps      = 100000000
nstxout-compressed = 50000
nstlog      = 10000
nstenergy   = 10000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
vdwtype     = Cut-off
rvdw        = 1.0
tcoupl      = V-rescale
tc-grps     = Protein Non-Protein
tau-t       = 0.1 0.1
ref-t       = 310 310
pcoupl      = no
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
"""

for mdp_name, mdp_content in [("em", mdp_em), ("nvt", mdp_nvt), ("npt", mdp_npt), ("md", mdp_md)]:
    for label in ["wt", "aib8"]:
        with open(f"{EXP_A}/gmx/{label}/{mdp_name}.mdp", "w") as f:
            f.write(mdp_content.strip() + "\n")

print("MDP files written")
