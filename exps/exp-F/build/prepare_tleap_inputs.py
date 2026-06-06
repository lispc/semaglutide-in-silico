#!/usr/bin/env python3
"""Prepare tleap inputs for exp-F: receptor+peptide PDB and aligned LNK mol2."""
import os

REPO = "/home/scroll/personal/semaglutide-in-silico"
BUILD = f"{REPO}/exps/exp-F/build"
STRUCT = f"{REPO}/exps/exp-F/structures"

combined_pdb = f"{STRUCT}/7ki0_receptor_plus_expD_ligand.pdb"
orig_mol2 = f"{REPO}/exps/exp-D/tleap/lnk_gglu_2oeg_pos.mol2"
out_pdb = f"{BUILD}/receptor_peptide.pdb"
out_mol2 = f"{BUILD}/lnk_aligned.mol2"

# 1. Write receptor_peptide.pdb (everything except LNK)
with open(combined_pdb) as f:
    lines = f.readlines()

with open(out_pdb, "w") as f:
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if "LNK" not in line[17:20]:
                f.write(line)
        elif line.startswith("TER"):
            f.write(line)
        elif line.startswith("END"):
            f.write(line)
    if not lines[-1].startswith("END"):
        f.write("END\n")

print(f"Wrote {out_pdb}")

# 2. Read LNK coordinates from combined PDB
pdb_coords = {}
with open(combined_pdb) as f:
    for line in f:
        if (line.startswith("ATOM") or line.startswith("HETATM")) and "LNK" in line[17:20]:
            name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            pdb_coords[name] = (x, y, z)

print(f"Read {len(pdb_coords)} LNK atoms from combined PDB")

# 3. Write aligned mol2: update coordinates for matching atoms, keep N from original
with open(orig_mol2) as f:
    mol2_lines = f.readlines()

updated_lines = []
in_atoms = False
atom_idx = 0
for line in mol2_lines:
    if "@<TRIPOS>ATOM" in line:
        in_atoms = True
        updated_lines.append(line)
        continue
    if "@<TRIPOS>BOND" in line or "@<TRIPOS>SUBSTRUCTURE" in line:
        in_atoms = False
        updated_lines.append(line)
        continue
    if in_atoms and line.strip():
        parts = line.split()
        name = parts[1]
        atom_idx += 1
        if name == "N":
            # Keep original N coordinate (will be removed by tleap)
            updated_lines.append(line)
        elif name in pdb_coords:
            x, y, z = pdb_coords[name]
            new_line = f"{atom_idx:>5} {name:<9} {x:9.4f} {y:9.4f} {z:9.4f} {parts[5]:>4} {parts[6]:>3} {parts[7]:<3} {parts[8]:>9}\n"
            updated_lines.append(new_line)
        else:
            print(f"Warning: atom {name} not found in PDB, keeping original")
            updated_lines.append(line)
    else:
        updated_lines.append(line)

with open(out_mol2, "w") as f:
    f.writelines(updated_lines)

print(f"Wrote {out_mol2}")

# 4. Verify alignment
print("\nVerification:")
with open(out_mol2) as f:
    lines = f.readlines()
    in_atoms = False
    count = 0
    for line in lines:
        if "@<TRIPOS>ATOM" in line:
            in_atoms = True
            continue
        if "@<TRIPOS>BOND" in line or "@<TRIPOS>SUBSTRUCTURE" in line:
            break
        if in_atoms and line.strip():
            count += 1
            if count <= 3:
                parts = line.split()
                print(f"  {parts[1]}: ({parts[2]}, {parts[3]}, {parts[4]})")
            elif count == 4:
                print("  ...")

print(f"Total atoms in aligned mol2: {count}")
