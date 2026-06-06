#!/usr/bin/env python3
"""Create a new mol2 file with N atom removed, NO BOND section."""
REPO = "/home/scroll/personal/semaglutide-in-silico"
MOL2_IN = f"{REPO}/exps/exp-F/build/lnk_aligned.mol2"
MOL2_OUT = f"{REPO}/exps/exp-F/build/lnk_stripped_nobond.mol2"

# Read atoms
atoms = []
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

# Remove N atom
new_atoms = [a for a in atoms if a['name'] != 'N']
print(f"After removing N: {len(new_atoms)} atoms")

# Write new mol2 without BOND section
with open(MOL2_OUT, "w") as f:
    f.write("@<TRIPOS>MOLECULE\n")
    f.write("LNK\n")
    f.write(f"{len(new_atoms):>5}     0     0     0     0\n")
    f.write("SMALL\n")
    f.write("GAFF2\n")
    f.write("linker-C18-noN\n")
    f.write("\n@<TRIPOS>ATOM\n")
    for i, a in enumerate(new_atoms):
        f.write(f"{i+1:>5} {a['name']:<9} {a['x']:9.4f} {a['y']:9.4f} {a['z']:9.4f} {a['type']:>4}   1 LNK {a['charge']:>9.6f}\n")
    f.write("\n@<TRIPOS>BOND\n")
    f.write("\n@<TRIPOS>SUBSTRUCTURE\n")
    f.write("     1 LNK         1 TEMP              0 ****  ****    0 ROOT\n")

print(f"Wrote {MOL2_OUT}")
