#!/usr/bin/env python3
"""Create a new mol2 file with N atom removed and explicit BOND section."""
import numpy as np

REPO = "/home/scroll/personal/semaglutide-in-silico"
MOL2_IN = f"{REPO}/exps/exp-F/build/lnk_aligned.mol2"
MOL2_OUT = f"{REPO}/exps/exp-F/build/lnk_stripped.mol2"

# Read atoms
atoms = []
bond_threshold = 1.8
with open(MOL2_IN) as f:
    lines = f.readlines()

in_atoms = False
for line in lines:
    if "@<TRIPOS>ATOM" in line:
        in_atoms = True
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
            'type': atype, 'charge': charge
        })

print(f"Read {len(atoms)} atoms")

# Find N atom to remove
n_idx = None
for i, a in enumerate(atoms):
    if a['name'] == 'N':
        n_idx = i
        break

if n_idx is None:
    raise ValueError("N atom not found")

print(f"Removing N atom at index {n_idx}")

# Compute bonds based on distance, filtering out H-H bonds
bonds = []
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        a1 = atoms[i]
        a2 = atoms[j]
        d = np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)
        if d < bond_threshold:
            # Skip H-H bonds
            is_h1 = a1['type'] in ('hc', 'hn', 'ho', 'h1', 'h2', 'h3', 'h4', 'h5', 'ha', 'hp')
            is_h2 = a2['type'] in ('hc', 'hn', 'ho', 'h1', 'h2', 'h3', 'h4', 'h5', 'ha', 'hp')
            if is_h1 and is_h2:
                continue
            bonds.append((i, j))

print(f"Detected {len(bonds)} bonds (excluding H-H)")

# Remove N atom and bonds involving N
new_atoms = [a for i, a in enumerate(atoms) if i != n_idx]
new_bonds = [(i, j) for i, j in bonds if i != n_idx and j != n_idx]

# Remap indices
idx_map = {}
new_idx = 1
for i, a in enumerate(atoms):
    if i != n_idx:
        idx_map[i] = new_idx
        new_idx += 1

# Remap bonds
remapped_bonds = []
for i, j in new_bonds:
    remapped_bonds.append((idx_map[i], idx_map[j]))

print(f"New atoms: {len(new_atoms)}, New bonds: {len(remapped_bonds)}")

# Write new mol2
with open(MOL2_OUT, "w") as f:
    f.write("@<TRIPOS>MOLECULE\n")
    f.write("LNK\n")
    f.write(f"{len(new_atoms):>5} {len(remapped_bonds):>5}     0     0     0\n")
    f.write("SMALL\n")
    f.write("GAFF2\n")
    f.write("linker-C18-noN\n")
    f.write("\n@<TRIPOS>ATOM\n")
    for i, a in enumerate(new_atoms):
        f.write(f"{i+1:>5} {a['name']:<9} {a['x']:9.4f} {a['y']:9.4f} {a['z']:9.4f} {a['type']:>4}   1 LNK {a['charge']:>9.6f}\n")
    f.write("\n@<TRIPOS>BOND\n")
    for i, (b1, b2) in enumerate(remapped_bonds):
        f.write(f"{i+1:>5} {b1:>5} {b2:>5} 1\n")
    f.write("\n@<TRIPOS>SUBSTRUCTURE\n")
    f.write("     1 LNK         1 TEMP              0 ****  ****    0 ROOT\n")

print(f"Wrote {MOL2_OUT}")
