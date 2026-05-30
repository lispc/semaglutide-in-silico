#!/usr/bin/env python3
"""Build LYA residue with proper AMBER backbone atom names, bypass prepgen."""
import os, subprocess, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

os.chdir("/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap")

# Step 1: Build ACE-LYA-NME
smiles = ("CC(=O)NC(CCCCNC(=O)CCC(N)C(=O)NCCOCCOCC(=O)NCCOCCOCC(=O)"
          "CCCCCCCCCCCCCCCCC(=O)[O-])C(=O)NC")
mol = Chem.MolFromSmiles(smiles); mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer(); n_atoms = mol.GetNumAtoms()
print(f"Built: {n_atoms} atoms")

# Step 2: Identify backbone atoms
# From SMILES structure: CC(=O)N[C@@H](CCCCN...)C(=O)NC
# ACE = atoms 0,1,2 (CH3, C=O, =O)
# Backbone N = atom 3, CA = atom 4
# Backbone C = find: C in LYA core bonded to =O and N(NME)
# NME = last N and CH3

# ACE atoms: first 3 heavy (0,1,2) + their H atoms (61,62,63 on atom 0)
ace_atoms = {0, 1, 2}
for i in range(n_atoms):
    atom = mol.GetAtomWithIdx(i)
    if atom.GetSymbol() == 'H':
        for n in atom.GetNeighbors():
            if n.GetIdx() in (0,): ace_atoms.add(i)

# Find NME: C=57-O=58 is backbone C=O, N=59-C=60 is NME
# NME: N=59, C=60, H on N59, H on C60
nme_atoms = {59, 60}
for i in (59, 60):
    atom = mol.GetAtomWithIdx(i)
    for n in atom.GetNeighbors():
        if n.GetSymbol() == 'H': nme_atoms.add(n.GetIdx())

# Backbone atoms (from SMILES analysis)
bb_n = 3; bb_ca = 4; bb_c = 57; bb_o = 58

# CB is atom 5 (first side chain C, attached to CA)
bb_cb = 5

print(f"Backbone: N={bb_n} CA={bb_ca} C={bb_c} O={bb_o} CB={bb_cb}")
print(f"ACE cap: {sorted(ace_atoms)}, NME cap: {sorted(nme_atoms)}")

# Step 3: Build LYA core atoms (everything except ACE/NME caps)
lya_atom_indices = [i for i in range(n_atoms) if i not in ace_atoms and i not in nme_atoms]
# Map old index → new index (0-based in LYA core)
old_to_new = {old: new for new, old in enumerate(lya_atom_indices)}

# Step 4: Assign AMBER backbone names, keep RDKit names for side chain
def get_lya_name(old_idx):
    if old_idx == bb_n: return 'N'
    if old_idx == bb_ca: return 'CA'
    if old_idx == bb_c: return 'C'
    if old_idx == bb_o: return 'O'
    if old_idx == bb_cb: return 'CB'
    # H atoms on backbone N, CA
    atom = mol.GetAtomWithIdx(old_idx)
    for n in atom.GetNeighbors():
        if n.GetIdx() == bb_n and atom.GetSymbol() == 'H': return 'H'
        if n.GetIdx() == bb_ca and atom.GetSymbol() == 'H': return 'HA'
        if n.GetIdx() == bb_cb and atom.GetSymbol() == 'H':
            # Need to distinguish HB2 and HB3
            # Just assign sequentially
            return 'HB2'  # will be fixed below

    # For other atoms, use element+number naming from RDKit
    sym = atom.GetSymbol()
    return f"{sym}{old_to_new[old_idx]+1}"

# Fix HB naming: find the two H on CB
cb_h_atoms = []
for old_idx in lya_atom_indices:
    atom = mol.GetAtomWithIdx(old_idx)
    if atom.GetSymbol() == 'H':
        for n in atom.GetNeighbors():
            if n.GetIdx() == bb_cb:
                cb_h_atoms.append(old_idx)

# Assign atom properties for LYA core
lya_atoms = []
for new_idx, old_idx in enumerate(lya_atom_indices):
    atom = mol.GetAtomWithIdx(old_idx)
    pos = conf.GetAtomPosition(old_idx)
    name = get_lya_name(old_idx)

    # Override for CB H atoms
    if old_idx in cb_h_atoms:
        h_num = cb_h_atoms.index(old_idx) + 2  # HB2, HB3
        name = f"HB{h_num}"

    # GAFF2 type
    sym = atom.GetSymbol()
    if sym == 'C': atype = 'c3'
    elif sym == 'O': atype = 'o'
    elif sym == 'N': atype = 'n'
    else: atype = 'hc'

    lya_atoms.append({
        'name': name, 'x': pos.x, 'y': pos.y, 'z': pos.z,
        'type': atype, 'charge': 0.0, 'old_idx': old_idx
    })

print(f"LYA core: {len(lya_atoms)} atoms")
print(f"Backbone atoms: {[(a['name'],a['old_idx']) for a in lya_atoms if a['name'] in ('N','H','CA','HA','C','O','CB','HB2','HB3')]}")

# Step 5: Build bond table for LYA core
lya_bonds = []
for bond in mol.GetBonds():
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    if i in lya_atom_indices and j in lya_atom_indices:
        lya_bonds.append((old_to_new[i]+1, old_to_new[j]+1, int(bond.GetBondTypeAsDouble())))

print(f"Bonds in LYA core: {len(lya_bonds)}")

# Step 6: Write mol2 for AM1-BCC
with open("lya_core.mol2", 'w') as f:
    n = len(lya_atoms)
    f.write("@<TRIPOS>MOLECULE\nLYA\n")
    f.write(f" {n} {len(lya_bonds)} 0 0 0\nSMALL\nGAFF2\nLYA residue\n\n")
    f.write("@<TRIPOS>ATOM\n")
    for i, a in enumerate(lya_atoms):
        f.write(f"{i+1:6d} {a['name']:6s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} {a['type']:4s} 1 LYA 0.0000\n")
    f.write("@<TRIPOS>BOND\n")
    for bi, (i, j, order) in enumerate(lya_bonds):
        f.write(f"{bi+1:6d} {i:6d} {j:6d} {order}\n")
    f.write("@<TRIPOS>SUBSTRUCTURE\n1 LYA 1 TEMP 0 **** **** 0 ROOT\n")

print("Saved: lya_core.mol2 (ready for AM1-BCC)")

# Step 7: Run AM1-BCC
print("\nRunning AM1-BCC...")
result = subprocess.run(
    "source /home/scroll/miniforge3/etc/profile.d/conda.sh && conda activate cgas-md && "
    "cd /home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap && "
    "rm -f ANTECHAMBER* ATOMTYPE* sqm.* 2>/dev/null && "
    "antechamber -i lya_core.mol2 -fi mol2 -o lya_bcc.mol2 -fo mol2 "
    "-c bcc -at gaff2 -nc 0 -rn LYA -pf yes 2>&1 | tail -3",
    shell=True, capture_output=True, text=True, timeout=600, executable='/bin/bash')
print(result.stdout.strip())
if result.returncode != 0:
    print("FAILED. See sqm.out")

print("\nDone. Next: check lya_bcc.mol2, then tleap sequence build.")
