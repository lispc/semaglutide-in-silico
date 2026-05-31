#!/usr/bin/env python3
"""Fix exp-D prmtop: add missing LNK bonds and NZ-C bond via ParmEd."""
import os, sys
import parmed as pmd
from parmed.topologyobjects import Bond
import numpy as np

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
BUILD = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def read_mol2_bonds(path):
    """Parse mol2 bond table, return list of (atom1_name, atom2_name) tuples."""
    bonds = []
    in_bond = False
    with open(path) as f:
        for line in f:
            ls = line.strip()
            if '@<TRIPOS>BOND' in ls: in_bond = True; continue
            if '@<TRIPOS>SUBSTRUCTURE' in ls: break
            if in_bond:
                parts = ls.split()
                if len(parts) >= 4:
                    bonds.append((int(parts[1]), int(parts[2])))
    return bonds

def fix_topology(vname):
    """Fix prmtop: add LNK bonds, add NZ-C bond, remove HZ."""
    prmtop_path = f"{TLEAP}/{vname}.prmtop"
    inpcrd_path = f"{TLEAP}/{vname}.inpcrd"
    mol2_path = f"{TLEAP}/lnk_{vname}_pos.mol2"

    if not os.path.exists(prmtop_path):
        print(f"SKIP {vname}: no prmtop")
        return
    if not os.path.exists(mol2_path):
        print(f"SKIP {vname}: no mol2")
        return

    print(f"\n=== Fixing {vname} ===")
    amber = pmd.load_file(prmtop_path, inpcrd_path)

    # Find Lys26 NZ and LNK atoms
    nz = None
    lnk_atoms = []  # in order of appearance in topology
    hz_atoms = []
    for a in amber.atoms:
        if a.residue.name == 'LNK':
            lnk_atoms.append(a)
        if a.name == 'NZ' and a.residue.name == 'LYS':
            nz = a
        if a.name.startswith('HZ') and a.residue.name == 'LYS':
            hz_atoms.append(a)

    if not nz:
        print("  ERROR: NZ not found")
        return
    print(f"  NZ: idx={nz.idx}, LNK atoms: {len(lnk_atoms)}, HZ atoms: {len(hz_atoms)}")

    # Remove HZ atoms completely
    # ParmEd doesn't easily delete atoms, so we'll just delete the bonds
    # Actually, ParmEd can delete atoms via indexing
    hz_indices = [a.idx for a in hz_atoms]
    if hz_indices:
        keep = [i for i in range(len(amber.atoms)) if i not in hz_indices]
        amber = amber[keep]
        # Re-find NZ and LNK atoms after deletion
        nz = None
        lnk_atoms = []
        for a in amber.atoms:
            if a.residue.name == 'LNK':
                lnk_atoms.append(a)
            if a.name == 'NZ' and a.residue.name == 'LYS':
                nz = a
        if not nz:
            print("  ERROR: NZ lost after HZ removal")
            return
        print(f"  After HZ removal: {len(amber.atoms)} atoms, LNK: {len(lnk_atoms)}")

    # Read mol2 bonds
    mol2_bonds = read_mol2_bonds(mol2_path)
    print(f"  Mol2 bonds: {len(mol2_bonds)}")

    # Map mol2 atom indices to ParmEd LNK atom indices
    # The LNK atoms in ParmEd are in order, matching mol2 order
    # mol2 atom 1 → lnk_atoms[0], mol2 atom 2 → lnk_atoms[1], etc.

    # Add LNK internal bonds
    bonds_added = 0
    for a1_mol2, a2_mol2 in mol2_bonds:
        if a1_mol2 <= len(lnk_atoms) and a2_mol2 <= len(lnk_atoms):
            a1 = lnk_atoms[a1_mol2 - 1]  # 0-based
            a2 = lnk_atoms[a2_mol2 - 1]
            amber.bonds.append(Bond(a1, a2))
            bonds_added += 1
    print(f"  Added {bonds_added} LNK internal bonds")

    # Add NZ-C bond (mol2 atom 2 = C = first heavy atom after N)
    lnk_c = lnk_atoms[0]  # After stripping N, atom 1 is C (the amide carbonyl)
    # Check: is it really C?
    if lnk_c.name != 'C':
        print(f"  WARNING: first LNK atom is {lnk_c.name}, expected C")
    amber.bonds.append(Bond(nz, lnk_c))
    print(f"  Added NZ-C bond: {nz.name}-{lnk_c.name}")

    # Save fixed prmtop
    out_prmtop = f"{TLEAP}/{vname}.prmtop"
    out_inpcrd = f"{TLEAP}/{vname}.inpcrd"
    amber.save(out_prmtop, format="amber", overwrite=True)
    amber.save(out_inpcrd, format="rst7", overwrite=True)
    print(f"  Saved: {out_prmtop} ({os.path.getsize(out_prmtop)/1e6:.1f}MB)")

    # Verify
    amber2 = pmd.load_file(out_prmtop, out_inpcrd)
    nz2 = next(a for a in amber2.atoms if a.name == 'NZ' and a.residue.name == 'LYS')
    lnk2_ids = {a.idx for a in amber2.atoms if a.residue.name == 'LNK'}
    lnk_bonds = [b for b in amber2.bonds if b.atom1.idx in lnk2_ids and b.atom2.idx in lnk2_ids]
    nz_lnk = [b for b in amber2.bonds if (b.atom1 is nz2 and b.atom2.idx in lnk2_ids) or (b.atom2 is nz2 and b.atom1.idx in lnk2_ids)]
    print(f"  Verify: LNK bonds={len(lnk_bonds)}, NZ-LNK={len(nz_lnk)}")

    return True

if __name__ == "__main__":
    os.chdir(TLEAP)
    for vname in VARIANTS:
        fix_topology(vname)
    print("\nDone!")
