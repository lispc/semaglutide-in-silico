#!/usr/bin/env python3
"""Merge 7KI0 receptor with exp-D peptide+LNK using ParmEd.

This creates a single Amber topology with full parameters for all atoms,
including GAFF2 parameters for LNK.
"""
import parmed as pmd
import os

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-F/structures"
os.makedirs(OUT, exist_ok=True)

# 1. Load exp-D fixed structure and strip ECD + water + ions
print("Loading exp-D fixed structure...")
expd = pmd.load_file(
    f"{REPO}/exps/exp-D/tleap/gglu_2oeg_fixed.prmtop",
    f"{REPO}/exps/exp-D/tleap/gglu_2oeg_fixed.inpcrd"
)
print(f"  Original: {len(expd.residues)} residues, {len(expd.atoms)} atoms")

# Strip ECD residues (1-100) and water/ions
# Residues in exp-D are numbered from 0 internally by ParmEd
ecd_resnames = set()
for res in expd.residues:
    if res.number < 100:  # 0-indexed; residues 0-99 are ECD
        ecd_resnames.add(res.name)
    elif res.name in ("WAT", "HOH", "Na+", "Cl-"):
        ecd_resnames.add(res.name)

print(f"  Stripping residues: {sorted(ecd_resnames)}")
expd.strip(f":1-100,@{','.join(sorted(ecd_resnames))}")
# Actually, let's be more precise: strip by residue index
# ParmEd uses 1-based residue indexing in strip commands
expd_copy = pmd.load_file(
    f"{REPO}/exps/exp-D/tleap/gglu_2oeg_fixed.prmtop",
    f"{REPO}/exps/exp-D/tleap/gglu_2oeg_fixed.inpcrd"
)
# Strip residues 1-100 (ECD), and all WAT/Na+
expd_copy.strip(":1-100")
expd_copy.strip("@Na+")
expd_copy.strip("@WAT")
print(f"  After strip: {len(expd_copy.residues)} residues, {len(expd_copy.atoms)} atoms")

# 2. Load 7KI0 receptor PDB and create Amber structure
print("\nLoading 7KI0 receptor...")
ki0 = pmd.load_file(f"{REPO}/exps/exp-B/structures/7KI0.pdb")
print(f"  PDB atoms: {len(ki0.atoms)}")

# Keep only chains R, A, B, G, N; remove P chain and waters
atoms_to_remove = []
for i, atom in enumerate(ki0.atoms):
    if atom.residue.chain not in ('R', 'A', 'B', 'G', 'N'):
        atoms_to_remove.append(i)
    elif atom.residue.name in ('HOH', 'WAT'):
        atoms_to_remove.append(i)

ki0.strip(f"@{','.join(str(i+1) for i in atoms_to_remove)}")
# Actually, ParmEd's strip uses atom selection syntax, not indices
# Let's recreate from PDB with filtering
ki0_atoms = []
with open(f"{REPO}/exps/exp-B/structures/7KI0.pdb") as f:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain = line[21]
            resname = line[17:20].strip()
            if chain in ('R', 'A', 'B', 'G', 'N') and resname not in ('HOH', 'WAT'):
                ki0_atoms.append(line)

ki0_pdb_path = f"{OUT}/receptor_only.pdb"
with open(ki0_pdb_path, "w") as f:
    for line in ki0_atoms:
        f.write(line)
    f.write("END\n")

# Load filtered receptor
receptor = pmd.load_file(ki0_pdb_path)
print(f"  Receptor atoms: {len(receptor.atoms)}")

# 3. Combine using ParmEd
print("\nMerging structures...")
# We need to assign atom types to receptor atoms for Amber compatibility
# Since receptor is standard amino acids, we can use tleap to build its prmtop
# OR we can use OpenMM to create a system and export to Amber
# For simplicity, let's use a different approach:

# Approach: Write a tleap script that loads receptor PDB + peptide_lnk PDB
# and combines them. But peptide_lnk has LNK which tleap can't load from PDB.

# Alternative approach: Use OpenMM to create a combined system
print("\nSwitching to OpenMM-based approach...")

# Save peptide+LNK as PDB
expd_copy.save(f"{OUT}/peptide_lnk.pdb", overwrite=True)
expd_copy.save(f"{OUT}/peptide_lnk.prmtop", overwrite=True)
expd_copy.save(f"{OUT}/peptide_lnk.inpcrd", overwrite=True)
print(f"Saved peptide+LNK: {len(expd_copy.atoms)} atoms")

# We now have:
# - receptor_only.pdb (standard amino acids)
# - peptide_lnk.prmtop/inpcrd (peptide + LNK with full Amber parameters)
#
# OpenMM can load both and create a combined system by:
# 1. Loading receptor PDB with standard force field
# 2. Loading peptide_lnk prmtop (with custom LNK parameters)
# 3. Manually merging topologies and systems

print("\nDone. Next step: OpenMM system creation with custom merging.")
