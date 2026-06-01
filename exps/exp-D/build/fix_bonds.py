#!/usr/bin/env python3
"""Fix exp-D prmtop: add missing LNK bonds + NZ-C bond using proper BondType."""
import os, sys
import parmed as pmd
from parmed.topologyobjects import Bond, BondType
import numpy as np

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
BUILD = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def read_mol2_bonds(path):
    """Parse positioned mol2, return list of (a1_name, a2_name, a1_idx, a2_idx)."""
    atoms = {}
    bonds = []
    in_atom = in_bond = False
    with open(path) as f:
        for line in f:
            ls = line.strip()
            if '@<TRIPOS>ATOM' in ls: in_atom = True; in_bond = False; continue
            if '@<TRIPOS>BOND' in ls: in_atom = False; in_bond = True; continue
            if in_atom and ls and not ls.startswith('@'):
                p = ls.split()
                if len(p) >= 6:
                    atoms[int(p[0])] = p[1]
            elif in_bond and ls and not ls.startswith('@'):
                p = ls.split()
                if len(p) >= 4 and p[0].isdigit():
                    bonds.append((int(p[1]), int(p[2])))
    return atoms, bonds

def fix_bonds(vname):
    """Fix prmtop: add LNK bonds + NZ-C bond."""
    prmtop = f"{TLEAP}/{vname}.prmtop"
    inpcrd = f"{TLEAP}/{vname}.inpcrd"
    if not os.path.exists(prmtop) or os.path.getsize(prmtop) == 0:
        print(f"SKIP {vname}: no valid prmtop")
        return

    print(f"\n=== {vname} ===")
    amber = pmd.load_file(prmtop, inpcrd)

    # Find NZ and HZ atoms
    nz = None; hz_idx = []
    for i, a in enumerate(amber.atoms):
        if a.name == 'NZ' and a.residue.name == 'LYS': nz = a
        if a.name.startswith('HZ') and a.residue.name == 'LYS': hz_idx.append(i)

    if not nz:
        print("  NZ not found!"); return

    # Find LNK atoms
    lnk_atoms = [a for a in amber.atoms if a.residue.name == 'LNK']
    print(f"  NZ idx={nz.idx}, HZ={len(hz_idx)}, LNK atoms={len(lnk_atoms)}")

    # Remove HZ atoms
    if hz_idx:
        keep = [i for i in range(len(amber.atoms)) if i not in hz_idx]
        amber = amber[keep]
        # Re-find NZ and LNK after deletion
        nz = next((a for a in amber.atoms if a.name == 'NZ' and a.residue.name == 'LYS'), None)
        lnk_atoms = [a for a in amber.atoms if a.residue.name == 'LNK']
        if not nz:
            print("  NZ lost!"); return
        print(f"  After HZ removal: {len(amber.atoms)} atoms")

    # Generate LNK bonds from distance-based detection
    # (positioned mol2 has no bond table)
    print("  Detecting LNK bonds from interatomic distances...")
    lnk_pos = np.array([[a.xx, a.xy, a.xz] for a in lnk_atoms])
    lnk_bonds = []
    # Bond length ranges per element pair (min, max)
    bond_ranges = {
        ('C','C'): (1.35, 1.65), ('C','O'): (1.15, 1.55),
        ('C','N'): (1.30, 1.55), ('C','H'): (0.95, 1.20),
        ('N','H'): (0.90, 1.15), ('O','H'): (0.85, 1.10),
    }
    for i in range(len(lnk_atoms)):
        for j in range(i+1, len(lnk_atoms)):
            e1 = lnk_atoms[i].element_name
            e2 = lnk_atoms[j].element_name
            key = tuple(sorted([e1, e2]))
            dmin, dmax = bond_ranges.get(key, bond_ranges.get((e1, e2), (0.5, 2.0)))
            d = np.linalg.norm(lnk_pos[i] - lnk_pos[j])
            if dmin < d < dmax:
                lnk_bonds.append((i+1, j+1))
    print(f"  Found {len(lnk_bonds)} distance-based bonds")

    # Bond types for different atom combinations
    # N3-c: amide (k=427, r=1.38)
    # c3-c3: C-C single (k=310, r=1.53)
    # c-o: carbonyl (k=570, r=1.23)
    # c3-n: C-N single (k=337, r=1.45)
    # c3-o: C-O single (k=320, r=1.41)
    bond_params = {
        ('N3', 'c'): (427.0, 1.38),
        ('c', 'N3'): (427.0, 1.38),
        ('c', 'c3'): (310.0, 1.53),
        ('c3', 'c'): (310.0, 1.53),
        ('c', 'o'): (570.0, 1.23),
        ('o', 'c'): (570.0, 1.23),
        ('c3', 'c3'): (310.0, 1.53),
        ('c3', 'n'): (337.0, 1.45),
        ('n', 'c3'): (337.0, 1.45),
        ('c3', 'o'): (320.0, 1.41),
        ('o', 'c3'): (320.0, 1.41),
        ('c3', 'hc'): (340.0, 1.09),
        ('hc', 'c3'): (340.0, 1.09),
        ('n', 'c'): (427.0, 1.38),
        ('c', 'n'): (427.0, 1.38),
        ('n', 'hn'): (434.0, 1.01),
        ('hn', 'n'): (434.0, 1.01),
    }

    bonds_added = 0
    for a1i, a2i in lnk_bonds:
        if a1i <= len(lnk_atoms) and a2i <= len(lnk_atoms):
            a1 = lnk_atoms[a1i - 1]
            a2 = lnk_atoms[a2i - 1]
            key = (a1.type, a2.type)
            k, r = bond_params.get(key, (300.0, 1.40))  # default
            bt = BondType(k, r)
            amber.bonds.append(Bond(a1, a2, type=bt))
            bonds_added += 1
    print(f"  Added {bonds_added} LNK internal bonds")

    # Add NZ-C bond (amide)
    lnk_c = lnk_atoms[0]  # First LNK atom = amide C
    bt_nzc = BondType(427.0, 1.38)
    amber.bonds.append(Bond(nz, lnk_c, type=bt_nzc))
    print(f"  Added NZ-C bond: {nz.name}-{lnk_c.name}")

    # Save
    out_prmtop = f"{TLEAP}/{vname}_fixed.prmtop"
    out_inpcrd = f"{TLEAP}/{vname}_fixed.inpcrd"
    try:
        amber.save(out_prmtop, format="amber", overwrite=True)
        amber.save(out_inpcrd, format="rst7", overwrite=True)
        print(f"  Saved: {out_prmtop} ({os.path.getsize(out_prmtop)/1e6:.1f}MB)")
    except Exception as e:
        print(f"  Save failed: {e}")

    # Also save as PDB for verification
    amber.save(f"{TLEAP}/{vname}_fixed.pdb", overwrite=True)

if __name__ == "__main__":
    os.chdir(TLEAP)
    for v in VARIANTS:
        fix_bonds(v)
    print("\nDone!")
