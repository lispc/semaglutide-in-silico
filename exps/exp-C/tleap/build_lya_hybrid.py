#!/usr/bin/env python3
"""Build LYA residue: ff14SB backbone + linker AM1-BCC side chain charges."""
import os, numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

os.chdir("/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap")

# Step 1: Build LYA core mol2 (same as before)
smiles = ("CC(=O)NC(CCCCNC(=O)CCC(N)C(=O)NCCOCCOCC(=O)NCCOCCOCC(=O)"
          "CCCCCCCCCCCCCCCCC(=O)[O-])C(=O)NC")
mol = Chem.MolFromSmiles(smiles); mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
AllChem.MMFFOptimizeMolecule(mol)
conf = mol.GetConformer(); n_atoms = mol.GetNumAtoms()

# ACE and NME caps
ace_atoms = {0, 1, 2}
for i in range(n_atoms):
    atom = mol.GetAtomWithIdx(i)
    if atom.GetSymbol() == 'H':
        for n in atom.GetNeighbors():
            if n.GetIdx() == 0: ace_atoms.add(i)

nme_atoms = {59, 60}
for i in (59, 60):
    atom = mol.GetAtomWithIdx(i)
    for n in atom.GetNeighbors():
        if n.GetSymbol() == 'H': nme_atoms.add(n.GetIdx())

# Backbone: N=3, CA=4, C=57, O=58, CB=5
bb_n=3; bb_ca=4; bb_c=57; bb_o=58; bb_cb=5

# LYA core indices
lya_indices = [i for i in range(n_atoms) if i not in ace_atoms and i not in nme_atoms]
old_to_new = {old: new for new, old in enumerate(lya_indices)}

# ff14SB Lys backbone charges (standard AMBER values)
FF14SB_CHARGES = {
    'N': -0.4157, 'H': 0.2719, 'CA': -0.0252, 'HA': 0.1079,
    'C': 0.5973, 'O': -0.5679,
    'CB': -0.0094, 'HB2': 0.0367, 'HB3': 0.0367,
}

# Load linker AM1-BCC charges (from charged form)
linker_charges = {}
with open("/tmp/linker_charged_bcc.mol2") as f:
    in_atom = False
    for line in f:
        if "@<TRIPOS>ATOM" in line: in_atom = True; continue
        if "@<TRIPOS>BOND" in line: break
        if in_atom and line.strip():
            p = line.split()
            if len(p) >= 9: linker_charges[p[1]] = float(p[-1])

print(f"Linker charges: {len(linker_charges)} atoms, total={sum(linker_charges.values()):.4f}")

# Build LYA atom list with charges
lya_atoms = []
for new_idx, old_idx in enumerate(lya_indices):
    atom = mol.GetAtomWithIdx(old_idx)
    pos = conf.GetAtomPosition(old_idx)
    sym = atom.GetSymbol()

    # Determine atom name
    if old_idx == bb_n: name = 'N'
    elif old_idx == bb_ca: name = 'CA'
    elif old_idx == bb_c: name = 'C'
    elif old_idx == bb_o: name = 'O'
    elif old_idx == bb_cb: name = 'CB'
    elif sym == 'H':
        # Find which heavy atom this H is attached to
        for n in atom.GetNeighbors():
            if n.GetIdx() == bb_n: name = 'H'; break
            elif n.GetIdx() == bb_ca: name = 'HA'; break
            elif n.GetIdx() == bb_cb:
                # HB2 or HB3: assign based on which H atom on CB
                cb_h_count = sum(1 for a in lya_atoms if a['name'].startswith('HB'))
                name = f"HB{cb_h_count+2}"
                break
        else:
            name = f"H{new_idx+1}"
    else:
        name = f"{sym}{new_idx+1}"

    # GAFF2 type
    if sym == 'C': atype = 'c3'
    elif sym == 'O': atype = 'o'
    elif sym == 'N': atype = 'n'
    else: atype = 'hc'

    # Charge: ff14SB for backbone, AM1-BCC for side chain
    if name in FF14SB_CHARGES:
        charge = FF14SB_CHARGES[name]
    elif old_idx <= 8:
        # Lys side chain CH2 groups (old_idx 5,6,7,8 and their H)
        # Use ff14SB Lys charges
        lys_sc_charges = {'CG': -0.0025, 'HG2': 0.0367, 'HG3': 0.0367,
                          'CD': -0.0025, 'HD2': 0.0367, 'HD3': 0.0367,
                          'CE': 0.1355, 'HE2': 0.0367, 'HE3': 0.0367,
                          'NZ': -0.3854}
        # Map based on position in side chain
        if old_idx == 5: charge = lys_sc_charges['CG']
        elif old_idx == 6: charge = lys_sc_charges['CD']
        elif old_idx == 7: charge = lys_sc_charges['CE']
        elif old_idx == 8: charge = lys_sc_charges['NZ']
        elif sym == 'H':
            for n in atom.GetNeighbors():
                if n.GetIdx() == 5: charge = 0.0367; break
                elif n.GetIdx() == 6: charge = 0.0367; break
                elif n.GetIdx() == 7: charge = 0.0367; break
                elif n.GetIdx() == 8: charge = 0.3000; break  # NZ H (will be removed)
            else: charge = 0.05
        else: charge = 0.05
    else:
        # Linker-FA atoms: map by approximate element position
        # Find matching atom in linker by element type and connectivity
        charge = 0.05  # fallback
        # Try to match by position in the chain
        for lname, lq in linker_charges.items():
            if lname[0] == sym and abs(new_idx - len(lya_indices) * 0.5) > 0:
                # Crude matching: use the average charge for this element type
                pass
        # Better: use a simple mapping
        if sym == 'C': charge = -0.10
        elif sym == 'O': charge = -0.55
        elif sym == 'N': charge = -0.50
        elif sym == 'H': charge = 0.05

    lya_atoms.append({
        'name': name, 'x': pos.x, 'y': pos.y, 'z': pos.z,
        'type': atype, 'charge': charge, 'old_idx': old_idx
    })

# Apply linker AM1-BCC charges to side chain atoms (old_idx >= 9)
# The side chain starts at old_idx 9 (the amide N connecting to γGlu)
# Map each side chain atom to the corresponding linker atom
# Linker atoms are numbered: O1, C2, O3, C4, C5, C6, N7, C8, ...
# Side chain atoms (old_idx 9+) should map to linker atoms

# Actually, the simplest approach: the linker portion is identical.
# Side chain starts at old_idx 9 = N (γGlu amide N)
# This corresponds to the N atom in the linker that's bonded to the amide C
# Let me just use the linker charges as-is for the side chain portion

# Count how many atoms are in the side chain (old_idx >= 9)
sc_atoms = [a for a in lya_atoms if a['old_idx'] >= 9]
linker_charge_list = list(linker_charges.values())
print(f"Side chain atoms: {len(sc_atoms)}, linker charges: {len(linker_charge_list)}")

# The side chain should have the same number of atoms as the linker
# If counts match, directly assign
if len(sc_atoms) == len(linker_charge_list):
    for i, a in enumerate(sc_atoms):
        a['charge'] = linker_charge_list[i]
    print("Directly assigned linker charges to side chain")
else:
    print(f"WARNING: size mismatch ({len(sc_atoms)} vs {len(linker_charge_list)})")
    # Adjust by adding/removing H atoms on NZ
    # The LYA has Lys NZ + linker, the free linker has γGlu N-term
    # The difference is 3 H atoms on NZ (which we'll remove for the amide bond)

total_q = sum(a['charge'] for a in lya_atoms)
print(f"Total LYA charge: {total_q:.4f}")

# Write mol2
with open("lya_hybrid.mol2", 'w') as f:
    n = len(lya_atoms)
    f.write("@<TRIPOS>MOLECULE\nLYA\n")
    f.write(f" {n} 0 0 0 0\nSMALL\nGAFF2\nLYA hybrid ff14SB+AM1-BCC\n\n")
    f.write("@<TRIPOS>ATOM\n")
    for i, a in enumerate(lya_atoms):
        f.write(f"{i+1:6d} {a['name']:6s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} {a['type']:4s} 1 LYA {a['charge']:10.4f}\n")

    # Bonds
    f.write("@<TRIPOS>BOND\n")
    bid = 0
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if i in lya_indices and j in lya_indices:
            bid += 1
            f.write(f"{bid:6d} {old_to_new[i]+1:6d} {old_to_new[j]+1:6d} {int(bond.GetBondTypeAsDouble())}\n")
    f.write("@<TRIPOS>SUBSTRUCTURE\n1 LYA 1 TEMP 0 **** **** 0 ROOT\n")

print(f"Saved: lya_hybrid.mol2 ({n} atoms)")
print("Next: test in tleap as standalone residue, then build full peptide sequence")
