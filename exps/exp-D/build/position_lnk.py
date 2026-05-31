#!/usr/bin/env python3
"""Position stripped LNK mol2 near Lys26 NZ, output positioned mol2.
Reads stripped mol2, applies Kabsch rotation+translation so amide C is near NZ.
"""
import numpy as np, os, shutil

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
BUILD = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

# Known Lys26 NZ position from ecd_pep.pdb
# From earlier: NZ at (-7.571, -0.448, 1.906)
NZ_POS = np.array([-7.571, -0.448, 1.906])
# CE at (-6.269, -0.196, 1.405)
CE_POS = np.array([-6.269, -0.196, 1.405])
# Direction CE→NZ
CE_TO_NZ = NZ_POS - CE_POS
CE_TO_NZ_DIR = CE_TO_NZ / np.linalg.norm(CE_TO_NZ)
print(f"CE→NZ: dir=({CE_TO_NZ_DIR[0]:.3f},{CE_TO_NZ_DIR[1]:.3f},{CE_TO_NZ_DIR[2]:.3f}), len={np.linalg.norm(CE_TO_NZ):.2f}")

def read_mol2(path):
    """Parse mol2 file, return atoms dict and bonds list."""
    atoms = {}
    bonds = []
    header = []
    in_atom = in_bond = False
    with open(path) as f:
        for line in f:
            line_s = line.strip()
            if '@<TRIPOS>ATOM' in line_s:
                in_atom = True; in_bond = False; header.append(line); continue
            if '@<TRIPOS>BOND' in line_s:
                in_atom = False; in_bond = True; header.append(line); continue
            if '@<TRIPOS>SUBSTRUCTURE' in line_s:
                in_atom = False; in_bond = False; header.append(line); continue
            parts = line_s.split()
            if in_atom and len(parts) >= 9:
                aid = int(parts[0])
                atoms[aid] = {'name': parts[1], 'x': float(parts[2]), 'y': float(parts[3]),
                              'z': float(parts[4]), 'atype': parts[5], 'resid': parts[6],
                              'resname': parts[7], 'charge': parts[8]}
            elif in_bond and len(parts) >= 4:
                bonds.append((int(parts[1]), int(parts[2]), parts[3]))
    return atoms, bonds

def write_mol2(path, atoms, bonds, natoms, nbonds):
    """Write positioned mol2."""
    # Renumber atoms sequentially
    old_to_new = {}
    for new_id, old_id in enumerate(sorted(atoms.keys()), 1):
        old_to_new[old_id] = new_id

    with open(path, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\nLNK\n")
        f.write(f" {natoms:5d} {nbonds:5d} 0 0 0\nSMALL\nGAFF2\npositioned linker\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for old_id in sorted(atoms.keys()):
            new_id = old_to_new[old_id]
            a = atoms[old_id]
            f.write(f"{new_id:6d} {a['name']:6s} {a['x']:10.4f} {a['y']:10.4f} "
                    f"{a['z']:10.4f} {a['atype']:4s} {a['resid']} {a['resname']} {a['charge']}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a1, a2, bt) in enumerate(bonds, 1):
            if a1 in old_to_new and a2 in old_to_new:
                f.write(f"{i:6d} {old_to_new[a1]:6d} {old_to_new[a2]:6d} {bt}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n1 LNK 1 TEMP 0 **** **** 0 ROOT\n")

def position_lnk(vname):
    """Position LNK mol2 near Lys26 NZ."""
    src = f"{BUILD}/{vname}_stripped.mol2"
    if not os.path.exists(src):
        print(f"SKIP {vname}: {src} not found")
        return None

    atoms, bonds = read_mol2(src)
    print(f"\n{vname}: {len(atoms)} atoms")

    # Find amide N (atom 1) and amide C (atom 2) in stripped mol2
    # The stripped mol2 starts with: N, C, O, then C13...
    keys = sorted(atoms.keys())
    n_id = keys[0]  # amide N
    c_id = keys[1]  # amide C
    o_id = keys[2]  # amide O

    # Get current N and C positions
    n_pos = np.array([atoms[n_id]['x'], atoms[n_id]['y'], atoms[n_id]['z']])
    c_pos = np.array([atoms[c_id]['x'], atoms[c_id]['y'], atoms[c_id]['z']])
    n_to_c = c_pos - n_pos
    print(f"  N→C: {np.linalg.norm(n_to_c):.2f} A")

    # Strategy:
    # Target is to place amide N→C along peptide CE→NZ direction
    # Place N at NZ position, extend C in CE→NZ direction (away from CE)
    # So C = NZ + CE_TO_NZ_DIR * |N→C|
    target_n = NZ_POS
    target_c = NZ_POS + CE_TO_NZ_DIR * np.linalg.norm(n_to_c)

    # Compute rotation from n_to_c to (target_c - target_n)
    v1 = n_to_c / np.linalg.norm(n_to_c)
    v2 = CE_TO_NZ_DIR
    cos_a = np.dot(v1, v2)
    if abs(cos_a - 1.0) < 1e-8:
        R = np.eye(3)
    elif abs(cos_a + 1.0) < 1e-8:
        R = -np.eye(3)
    else:
        k = np.cross(v1, v2)
        k = k / np.linalg.norm(k)
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R = np.eye(3) + np.sin(np.arccos(cos_a)) * K + (1 - cos_a) * K @ K

    # Apply: rotate around N, then translate N to NZ
    for aid in atoms:
        pos = np.array([atoms[aid]['x'], atoms[aid]['y'], atoms[aid]['z']])
        rel = pos - n_pos
        new_pos = R @ rel + target_n
        atoms[aid]['x'], atoms[aid]['y'], atoms[aid]['z'] = new_pos

    # Verify
    new_n = np.array([atoms[n_id]['x'], atoms[n_id]['y'], atoms[n_id]['z']])
    new_c = np.array([atoms[c_id]['x'], atoms[c_id]['y'], atoms[c_id]['z']])
    nz_to_newc = new_c - NZ_POS
    print(f"  NZ→C distance: {np.linalg.norm(nz_to_newc):.2f} A")
    print(f"  NZ→C dir vs CE→NZ: dot={np.dot(nz_to_newc/np.linalg.norm(nz_to_newc), CE_TO_NZ_DIR):.3f}")

    # Write positioned mol2
    out_path = f"{TLEAP}/lnk_{vname}_pos.mol2"
    write_mol2(out_path, atoms, bonds, len(atoms), len(bonds))
    print(f"  Saved: {out_path}")
    return out_path

if __name__ == "__main__":
    os.makedirs(TLEAP, exist_ok=True)
    for v in VARIANTS:
        position_lnk(v)
    print("\nDone!")
