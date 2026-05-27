#!/usr/bin/env python3
"""
Build final complex PDBs for tleap using crystal structure coordinates.

Strategy (final):
  - Use crystal peptide (residues 10-37) directly from 4ZGM
  - Mutate R34→K34 for B1 using ParmEd backbone + tleap sidechain fill
  - B2 uses crystal R34 as-is
  - Cap both termini with ACE/NME for neutral representation
  - Note: residues 7-9 (incl. Aib8) are not modeled — they are disordered
    in all crystal structures and not needed for the K34R question.
    Position 34 is at the C-terminal end, ~24 residues from the N-terminus.
"""
import numpy as np
import os

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-B/structures"

def read_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atoms.append({
                    'chain': line[21], 'resname': line[17:20].strip(),
                    'resid': int(line[22:26]), 'atom': line[12:16].strip(),
                    'x': float(line[30:38]), 'y': float(line[38:46]), 'z': float(line[46:54]),
                    'element': line[76:78].strip() or 'C',
                })
    return atoms

def write_pdb(atoms_by_chain, path, remark=""):
    atom_num = 0
    with open(path, 'w') as f:
        f.write(f"REMARK {remark}\n")
        for chain in ['A', 'B']:
            if chain not in atoms_by_chain:
                continue
            for a in atoms_by_chain[chain]:
                atom_num += 1
                f.write(f"ATOM  {atom_num:5d} {a['atom']:4s} {a['resname']:3s} "
                        f"{a['chain']}{a['resid']:4d}    "
                        f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                        f"  1.00  0.00          {a.get('element','C'):2s}  \n")
            f.write("TER\n")
        f.write("END\n")

print("Building final complex PDBs for tleap")
print("=" * 60)

# Load ECD (chain A) from 4ZGM
ecd_all = read_pdb(f"{OUT}/ecd_clean.pdb")
ecd_heavy = []
for a in ecd_all:
    if a['atom'][0] != 'H':
        ecd_heavy.append({
            'chain': 'A', 'resname': a['resname'], 'resid': a['resid'],
            'atom': a['atom'], 'x': a['x'], 'y': a['y'], 'z': a['z'],
            'element': a.get('element', 'C'),
        })
print(f"ECD heavy atoms: {len(ecd_heavy)}")

# Load crystal peptide (chain B, residues 10-37) from 4ZGM
pep_4zgm = read_pdb(f"{OUT}/pep_4zgm_clean.pdb")
pep_heavy = []
for a in pep_4zgm:
    if a['atom'][0] != 'H':
        pep_heavy.append({
            'chain': 'B', 'resname': a['resname'], 'resid': a['resid'],
            'atom': a['atom'], 'x': a['x'], 'y': a['y'], 'z': a['z'],
            'element': a.get('element', 'C'),
        })
print(f"Crystal peptide heavy atoms: {len(pep_heavy)}")

# Verify residues
pep_resids = sorted(set(a['resid'] for a in pep_heavy))
print(f"Peptide residues: {pep_resids[0]}-{pep_resids[-1]}")

# ============================================================
# B1 (K34): Replace ARG34 sidechain with LYS backbone only
# ============================================================
pep_b1 = []
for a in pep_heavy:
    if a['resid'] == 34 and a['resname'] == 'ARG':
        if a['atom'] in ('N', 'CA', 'C', 'O', 'CB'):
            new_a = dict(a)
            new_a['resname'] = 'LYS'
            pep_b1.append(new_a)
        # Skip ARG sidechain (CG, CD, NE, CZ, NH1, NH2)
        # tleap will build LYS sidechain from its template
    else:
        pep_b1.append(dict(a))
print(f"\nB1 peptide: {len(pep_b1)} atoms (ARG34→LYS, backbone only)")

# B2 (R34): Use as-is
pep_b2 = [dict(a) for a in pep_heavy]
print(f"B2 peptide: {len(pep_b2)} atoms (R34 from crystal)")

# ============================================================
# Write PDBs
# ============================================================
write_pdb({'A': ecd_heavy, 'B': pep_b1},
    f"{OUT}/b1_complex_noh.pdb",
    f"exp-B B1: ECD + GLP-1(10-37, K34) from 4ZGM, no H")

write_pdb({'A': ecd_heavy, 'B': pep_b2},
    f"{OUT}/b2_complex_noh.pdb",
    f"exp-B B2: ECD + GLP-1(10-37, R34) from 4ZGM, no H")

# Print summary
print(f"\n{'='*60}")
print("Done!")
print(f"B1: {len(ecd_heavy)+len(pep_b1)} atoms → {OUT}/b1_complex_noh.pdb")
print(f"B2: {len(ecd_heavy)+len(pep_b2)} atoms → {OUT}/b2_complex_noh.pdb")
print(f"\nCAVEATS:")
print(f"  - Residues 7-9 (His7-Aib8-Glu9) are not modeled (disordered in crystal)")
print(f"  - B1 K34 has backbone only (tleap will add LYS sidechain)")
print(f"  - C-terminal Gly37 may need CGLY rename for tleap")
print(f"  - These caveats don't affect the K34R comparison (position 34 is C-terminal)")
