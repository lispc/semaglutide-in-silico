#!/usr/bin/env python3
"""Strip LYA mol2: parse mol2 as text, keep only linker+tail atoms.
Output: C(=O)-linker-C18-COO⁻ for bonding to peptide Lys26 NZ.
"""
import os

BUILD = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"
TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"

# Heavy atoms to DELETE (and their H atoms):
# ACE: C1,C2,O3; Backbone: N4,C5; Sidechain: C6,C7,C8,C9
# NME: C32,O33,N34,C35
# KEEP: N10, C11, O12, + C13-C31 (linker + C18 chain + COO⁻)
DELETE_HEAVY = {'C1','C2','O3','N4','C5','C6','C7','C8','C9',
                'C32','O33','N34','C35'}
# N10→N, C11→C, O12→O
RENAME = {'N10': 'N', 'C11': 'C', 'O12': 'O'}

def gaff2_type(name, elem):
    """Assign GAFF2 atom type based on atom name/role in stripped LNK."""
    if name == 'N': return 'n'      # amide N (sp2)
    if name == 'C': return 'c'      # amide carbonyl C (sp2)
    if name == 'O': return 'o'      # amide carbonyl O
    # C13-C28: CH2 chain carbons → c3
    # C29: next to COO⁻ → c3
    # C30 (COO⁻ C): named C31 in variants? Need to detect.
    # Actually names vary by variant. Use element-based detection.
    if elem == 'C': return 'c3'     # all other C are sp3 CH2
    if elem.startswith('O'): return 'o'  # all O are carbonyl/carboxyl
    if elem.startswith('N'): return 'n'
    if elem.startswith('H'): return 'hc'
    return 'du'

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def strip_mol2(vname):
    """Parse mol2 text, strip unwanted atoms, renumber."""
    mol2_in = f"{BUILD}/{vname}_bcc.mol2"
    if not os.path.exists(mol2_in):
        return None

    with open(mol2_in) as f:
        lines = f.readlines()

    # Parse atoms from bcc mol2, bonds from raw mol2
    atoms = {}  # id -> {name, x, y, z, atype, resid, resname, charge}
    bonds = []  # [(a1, a2, btype)]
    in_atom = in_bond = False

    for line in lines:
        line = line.strip()
        if '@<TRIPOS>ATOM' in line:
            in_atom = True; in_bond = False; continue
        if '@<TRIPOS>BOND' in line:
            in_atom = False; in_bond = True; continue
        if '@<TRIPOS>SUBSTRUCTURE' in line:
            in_atom = False; in_bond = False; continue

        parts = line.split()
        if in_atom and len(parts) >= 9:
            aid = int(parts[0])
            atoms[aid] = {
                'name': parts[1], 'x': parts[2], 'y': parts[3], 'z': parts[4],
                'atype': parts[5], 'resid': parts[6], 'resname': parts[7],
                'charge': parts[8]
            }
        elif in_bond and len(parts) >= 4:
            bonds.append((int(parts[1]), int(parts[2]), parts[3]))

    # If bcc mol2 has no bonds, read bonds from raw mol2
    if len(bonds) == 0:
        raw_mol2 = f"{BUILD}/{vname}_raw.mol2"
        if os.path.exists(raw_mol2):
            with open(raw_mol2) as rf:
                in_bond_raw = False
                for rline in rf:
                    rs = rline.strip()
                    if '@<TRIPOS>BOND' in rs: in_bond_raw = True; continue
                    if '@<TRIPOS>SUBSTRUCTURE' in rs: break
                    if in_bond_raw:
                        rparts = rs.split()
                        if len(rparts) >= 4:
                            bonds.append((int(rparts[1]), int(rparts[2]), rparts[3]))
            print(f"  Loaded {len(bonds)} bonds from raw mol2")

    # Find which atom IDs to delete (heavy + their H atoms)
    delete_ids = set()
    for aid, a in atoms.items():
        if a['name'] in DELETE_HEAVY:
            delete_ids.add(aid)

    # Find H atoms bonded to deleted heavy atoms
    for a1, a2, _ in bonds:
        if a1 in delete_ids and a2 not in delete_ids:
            if atoms[a2]['name'].startswith('H'):
                delete_ids.add(a2)
        if a2 in delete_ids and a1 not in delete_ids:
            if atoms[a1]['name'].startswith('H'):
                delete_ids.add(a1)

    keep_ids = sorted(set(atoms.keys()) - delete_ids)

    # Build old→new atom ID mapping
    old_to_new = {}
    for new_id, old_id in enumerate(keep_ids, 1):
        old_to_new[old_id] = new_id

    # Rename atoms and fix types
    for old_id in keep_ids:
        old_name = atoms[old_id]['name']
        if old_name in RENAME:
            atoms[old_id]['name'] = RENAME[old_name]
        atoms[old_id]['atype'] = gaff2_type(atoms[old_id]['name'], atoms[old_id]['name'][0])

    # Renumber bonds (only bonds where both atoms are kept)
    new_bonds = []
    for a1, a2, bt in bonds:
        if a1 in old_to_new and a2 in old_to_new:
            new_bonds.append((old_to_new[a1], old_to_new[a2], bt))

    # Count total charge
    total_charge = sum(float(atoms[aid]['charge']) for aid in keep_ids)

    # Write new mol2
    out_mol2 = f"{BUILD}/{vname}_stripped.mol2"
    with open(out_mol2, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\nLNK\n")
        f.write(f" {len(keep_ids):5d} {len(new_bonds):5d} 0 0 0\n")
        f.write("SMALL\nGAFF2\nlinker-C18\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for new_id, old_id in enumerate(keep_ids, 1):
            a = atoms[old_id]
            f.write(f"{new_id:6d} {a['name']:6s} {float(a['x']):10.4f} "
                    f"{float(a['y']):10.4f} {float(a['z']):10.4f} "
                    f"{a['atype']:4s} {a['resid']} LNK {a['charge']}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a1, a2, bt) in enumerate(new_bonds, 1):
            f.write(f"{i:6d} {a1:6d} {a2:6d} {bt}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n1 LNK 1 TEMP 0 **** **** 0 ROOT\n")

    print(f"{vname}: {len(atoms)} → {len(keep_ids)} atoms (charge={total_charge:.3f})")
    return out_mol2, len(keep_ids)

if __name__ == "__main__":
    for vname in VARIANTS:
        result = strip_mol2(vname)
    print("Done!")
