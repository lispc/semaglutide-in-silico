#!/usr/bin/env python3
"""Fix hydrogen atom types in mol2: N-attached H should be 'hn', not 'hc'."""
import numpy as np

REPO = "/home/scroll/personal/semaglutide-in-silico"
MOL2_IN = f"{REPO}/exps/exp-F/build/lnk_stripped.mol2"
MOL2_OUT = f"{REPO}/exps/exp-F/build/lnk_stripped_fixed.mol2"

# Read atoms
atoms = []
with open(MOL2_IN) as f:
    lines = f.readlines()

in_atoms = False
atom_lines = []
for i, line in enumerate(lines):
    if "@<TRIPOS>ATOM" in line:
        in_atoms = True
        atom_lines.append(i)
        continue
    if "@<TRIPOS>BOND" in line or "@<TRIPOS>SUBSTRUCTURE" in line:
        in_atoms = False
        continue
    if in_atoms and line.strip():
        parts = line.split()
        idx = int(parts[0])
        name = parts[1]
        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
        atype = parts[5]
        charge = float(parts[8])
        atoms.append({
            'idx': idx, 'name': name, 'x': x, 'y': y, 'z': z,
            'type': atype, 'charge': charge, 'line_idx': i
        })

print(f"Read {len(atoms)} atoms")

# Read bonds
bonds = []
in_bonds = False
for line in lines:
    if "@<TRIPOS>BOND" in line:
        in_bonds = True
        continue
    if "@<TRIPOS>SUBSTRUCTURE" in line:
        break
    if in_bonds and line.strip():
        parts = line.split()
        b1, b2 = int(parts[1]), int(parts[2])
        bonds.append((b1, b2))

print(f"Read {len(bonds)} bonds")

# Find N-attached hydrogens and fix their type
fixed_count = 0
for b1, b2 in bonds:
    a1 = atoms[b1-1]  # 1-indexed to 0-indexed
    a2 = atoms[b2-1]
    # Check if one is nitrogen (type 'n') and the other is hydrogen (type 'hc')
    if a1['type'] == 'n' and a2['type'] == 'hc':
        print(f"Fixing {a2['name']} (attached to {a1['name']}): hc -> hn")
        a2['type'] = 'hn'
        fixed_count += 1
    elif a2['type'] == 'n' and a1['type'] == 'hc':
        print(f"Fixing {a1['name']} (attached to {a2['name']}): hc -> hn")
        a1['type'] = 'hn'
        fixed_count += 1

print(f"Fixed {fixed_count} hydrogen atom types")

# Write fixed mol2
with open(MOL2_OUT, "w") as f:
    in_atoms = False
    in_bonds = False
    for i, line in enumerate(lines):
        if "@<TRIPOS>ATOM" in line:
            in_atoms = True
            f.write(line)
            continue
        if "@<TRIPOS>BOND" in line:
            in_atoms = False
            in_bonds = True
            f.write(line)
            continue
        if "@<TRIPOS>SUBSTRUCTURE" in line:
            in_bonds = False
            f.write(line)
            continue
        if in_atoms and line.strip():
            # Find corresponding atom
            parts = line.split()
            idx = int(parts[0])
            atom = atoms[idx-1]
            new_line = f"{idx:>5} {atom['name']:<9} {atom['x']:9.4f} {atom['y']:9.4f} {atom['z']:9.4f} {atom['type']:>4}   1 LNK {atom['charge']:>9.6f}\n"
            f.write(new_line)
        else:
            f.write(line)

print(f"Wrote {MOL2_OUT}")
