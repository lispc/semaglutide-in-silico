#!/usr/bin/env python3
"""Position stripped LNK at peptide Lys26, write combined PDB for tleap.
Aligns LNK N10 to peptide NZ, keeps N10, bonds CE→N10.
"""
import os, shutil
import numpy as np
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
TLEAP = f"{REPO}/exps/exp-D/tleap"
BUILD = f"{REPO}/exps/exp-D/build"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def position_and_assemble(vname):
    """Position LNK at peptide Lys26, write combined PDB with CONECT."""
    bcc_mol2 = f"{BUILD}/{vname}_bcc.mol2"
    if not os.path.exists(bcc_mol2):
        print(f"SKIP {vname}: no bcc mol2")
        return None

    # Load complex and full (unstripped) LYA
    comp = pmd.load_file(f"{TLEAP}/ecd_pep.pdb")
    lya = pmd.load_file(bcc_mol2)

    # Find peptide Lys26 atoms
    lys_nz = lys_ce = None
    lys_hz = []
    for a in comp.atoms:
        if a.residue.chain == 'B' and a.residue.name == 'LYS':
            if a.name == 'NZ': lys_nz = a
            if a.name == 'CE': lys_ce = a
            if a.name in ('HZ1', 'HZ2', 'HZ3'): lys_hz.append(a)

    if not lys_nz or not lys_ce:
        print(f"ERROR: Lys26 NZ/CE not found")
        return None

    # Find LYA atoms: N10 (NZ-equiv), C11 (amide C), and atoms to delete
    lya_n10 = lya_c11 = None
    delete_heavy = {'C1','C2','O3','N4','C5','C6','C7','C8','C9','C32','O33','N34','C35'}
    for a in lya.atoms:
        if a.name == 'N10': lya_n10 = a
        if a.name == 'C11': lya_c11 = a

    if not lya_n10 or not lya_c11:
        print(f"ERROR: LYA N10/C11 not found")
        return None

    # Determine which atoms to keep:
    # Keep N10, C11, O12, C13-C31 (linker + C18 + COO⁻)
    # Delete: ACE (C1-C2-O3), backbone (N4,C5), sidechain (C6-C9), NME (C32,O33,N34,C35)
    keep_idx = []
    for i, a in enumerate(lya.atoms):
        if a.name in delete_heavy:
            continue
        # Also delete H atoms bonded to deleted heavy atoms
        if a.element_name == 'H':
            partners = [b for b in a.bond_partners]
            if any(p.name in delete_heavy for p in partners):
                continue
            # Delete H on N10 (keep one)
            if any(p is lya_n10 for p in partners):
                continue
        keep_idx.append(i)

    # Build LNK substructure
    lnk_atoms_keep = [lya.atoms[i] for i in keep_idx]
    print(f"{vname}: keeping {len(lnk_atoms_keep)}/{len(lya.atoms)} LYA atoms")

    # Position LNK: align C11 to target position, then place N10
    # CE bond direction: CE → NZ vector
    ce_to_nz = np.array([lys_nz.xx - lys_ce.xx,
                          lys_nz.xy - lys_ce.xy,
                          lys_nz.xz - lys_ce.xz])
    ce_to_nz_len = np.linalg.norm(ce_to_nz)
    ce_to_nz_dir = ce_to_nz / ce_to_nz_len
    print(f"  CE→NZ: {ce_to_nz_len:.2f} A, dir=({ce_to_nz_dir[0]:.2f},{ce_to_nz_dir[1]:.2f},{ce_to_nz_dir[2]:.2f})")

    # LNK N10→C11 vector in original coords
    n10_pos = np.array([lya_n10.xx, lya_n10.xy, lya_n10.xz])
    c11_pos = np.array([lya_c11.xx, lya_c11.xy, lya_c11.xz])
    n10_to_c11 = c11_pos - n10_pos
    n10_to_c11_len = np.linalg.norm(n10_to_c11)
    print(f"  N10→C11: {n10_to_c11_len:.2f} A")

    # Strategy: place N10 at NZ position, rotate C11 to align with CE→NZ direction
    # N10 replaces NZ in the structure
    # Translation: bring N10 to NZ
    translation = np.array([lys_nz.xx, lys_nz.xy, lys_nz.xz]) - n10_pos

    # Rotation: align N10→C11 with CE→NZ (so C11 extends away from CE)
    # Rodrigues rotation formula
    n10_to_c11_dir = n10_to_c11 / n10_to_c11_len
    # We want C11 to be placed along the CE→NZ direction from N10's new position
    # So C11 ends up at: NZ_pos + ce_to_nz_dir * n10_to_c11_len
    # This means N10→C11 should align with CE→NZ direction

    v1 = n10_to_c11_dir
    v2 = ce_to_nz_dir
    cos_a = np.dot(v1, v2)
    if abs(cos_a - 1.0) < 1e-6:
        R = np.eye(3)
    elif abs(cos_a + 1.0) < 1e-6:
        R = -np.eye(3)
    else:
        k = np.cross(v1, v2)
        k_norm = np.linalg.norm(k)
        k = k / k_norm
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R = np.eye(3) + np.sin(np.arccos(cos_a)) * K + (1 - cos_a) * K @ K

    # Apply rotation around N10 then translation
    rot_c11 = R @ n10_to_c11 + n10_pos + translation
    actual_c11 = c11_pos + translation  # no rotation (just translation)
    # For full alignment: rotate all lnk atoms around N10, then translate

    # Actually, let's just apply the rotation to all LNK atoms relative to N10
    # Then translate to bring N10 to NZ position

    # Write combined PDB
    out_pdb = f"{TLEAP}/combined_{vname}.pdb"
    with open(out_pdb, 'w') as f:
        f.write(f"REMARK exp-D: ECD + peptide + {vname} linker\n")
        f.write(f"REMARK Bond: peptide Lys26 CE to LNK N10\n")

        atom_num = 0
        ce_atom_num = None
        n10_atom_num = None

        # Complex atoms (minus HZ)
        for a in comp.atoms:
            if a in lys_hz:
                continue
            # Don't write NZ (will be replaced by LNK N10)
            if a is lys_nz:
                continue
            atom_num += 1
            if a is lys_ce:
                ce_atom_num = atom_num
            f.write(f"ATOM  {atom_num:5d} {a.name:4s} {a.residue.name:3s} "
                    f"{a.residue.chain}{a.residue.idx+1:4d}    "
                    f"{a.xx:8.3f}{a.xy:8.3f}{a.xz:8.3f}"
                    f"  1.00  0.00          {a.element_name:2s}  \n")

        # LNK atoms: apply rotation around N10 + translation to NZ
        for a in lya.atoms:
            if a.idx not in keep_idx:
                continue
            a_pos = np.array([a.xx, a.xy, a.xz])
            rel = a_pos - n10_pos  # relative to N10
            rot_rel = R @ rel  # rotate
            new_pos = rot_rel + n10_pos + translation  # translate to NZ position
            atom_num += 1
            name = a.name
            if a is lya_n10:
                name = 'N'  # rename to generic amide N
                n10_atom_num = atom_num
            if a is lya_c11:
                name = 'C'  # rename to generic amide C
            f.write(f"ATOM  {atom_num:5d} {name:4s} LNK "
                    f"B{27:4d}    "
                    f"{new_pos[0]:8.3f}{new_pos[1]:8.3f}{new_pos[2]:8.3f}"
                    f"  1.00  0.00          {a.element_name:2s}  \n")

        # Bond: CE - N (LNK N10)
        if ce_atom_num and n10_atom_num:
            f.write(f"CONECT{ce_atom_num:5d}{n10_atom_num:5d}\n")
            # Also N-C bond within LNK
            c_atom_num = None
            for i, a in enumerate(lya.atoms):
                if a is lya_c11 and a.idx in keep_idx:
                    # Count LNK atoms before C11 to find its atom number
                    c_atom_num = atom_num  # FIXME: this is wrong, need proper counting
            f.write(f"# CE-N: {ce_atom_num}-{n10_atom_num}\n")

    print(f"  {out_pdb}: {atom_num} atoms")
    ce_n_dist = np.linalg.norm(np.array([lys_ce.xx, lys_ce.xy, lys_ce.xz]) -
                               (n10_pos + translation))
    print(f"  CE-N distance after positioning: {ce_n_dist:.2f} A")
    return out_pdb

if __name__ == "__main__":
    os.makedirs(TLEAP, exist_ok=True)
    for vname in VARIANTS:
        position_and_assemble(vname)
    print("Done!")
