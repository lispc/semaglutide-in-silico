#!/usr/bin/env python3
"""Build combined ECD-peptide-LNK PDB with CONECT records for tleap.
tleap loads PDB + CONECT to establish the NZ-C bond.
"""
import os, shutil, subprocess
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
TLEAP = f"{REPO}/exps/exp-D/tleap"
BUILD = f"{REPO}/exps/exp-D/build"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def build_combined_pdb(vname):
    """Write combined PDB: ECD-peptide + LNK, with NZ-C CONECT bond."""
    stripped_mol2 = f"{BUILD}/{vname}_stripped.mol2"
    if not os.path.exists(stripped_mol2):
        print(f"SKIP {vname}: no stripped mol2")
        return None

    # Load complex and LNK
    complex_pdb = f"{TLEAP}/ecd_pep.pdb"
    comp = pmd.load_file(complex_pdb)
    lnk = pmd.load_file(stripped_mol2)

    # Find Lys26 and its NZ/HZ atoms
    # After ParmEd renumbering, residue numbers changed. Find by chain B, name LYS
    lys26_hz = []
    lys26_nz = None
    for a in comp.atoms:
        if a.residue.chain == 'B' and a.residue.name == 'LYS':
            if a.name in ('HZ1', 'HZ2', 'HZ3'):
                lys26_hz.append(a)
            if a.name == 'NZ':
                lys26_nz = a

    if not lys26_nz:
        print(f"ERROR: Lys26 NZ not found")
        return None

    # LNK amide C (renamed from C11 to 'C' in strip)
    lnk_c = None
    for a in lnk.atoms:
        if a.name == 'C':
            lnk_c = a
            break
    if not lnk_c:
        print(f"ERROR: LNK amide C not found")
        return None

    print(f"{vname}: NZ idx={lys26_nz.idx}, LNK C idx={lnk_c.idx}, "
          f"HZ to remove: {len(lys26_hz)}")
    print(f"  NZ pos: ({lys26_nz.xx:.1f},{lys26_nz.xy:.1f},{lys26_nz.xz:.1f})")
    print(f"  C pos:  ({lnk_c.xx:.1f},{lnk_c.xy:.1f},{lnk_c.xz:.1f})")
    dist = ((lys26_nz.xx-lnk_c.xx)**2 + (lys26_nz.xy-lnk_c.xy)**2 + (lys26_nz.xz-lnk_c.xz)**2)**0.5
    print(f"  NZ-C distance: {dist:.2f} A")

    # Build combined PDB
    out_pdb = f"{TLEAP}/combined_{vname}.pdb"
    with open(out_pdb, 'w') as f:
        f.write(f"REMARK exp-D: ECD + peptide + {vname} linker\n")
        f.write(f"REMARK Lys26 NZ bonded to LNK amide C via CONECT\n")

        atom_num = 0
        nz_atom_num = None
        lnk_c_atom_num = None

        # Write complex atoms (minus HZ)
        for a in comp.atoms:
            if a in lys26_hz:
                continue  # skip HZ atoms
            atom_num += 1
            if a is lys26_nz:
                nz_atom_num = atom_num
            f.write(f"ATOM  {atom_num:5d} {a.name:4s} {a.residue.name:3s} "
                    f"{a.residue.chain}{a.residue.idx+1:4d}    "
                    f"{a.xx:8.3f}{a.xy:8.3f}{a.xz:8.3f}"
                    f"  1.00  0.00          {a.element_name:2s}  \n")

        # Write LNK atoms
        for a in lnk.atoms:
            atom_num += 1
            if a is lnk_c:
                lnk_c_atom_num = atom_num
            f.write(f"ATOM  {atom_num:5d} {a.name:4s} LNK "
                    f"B{27:4d}    "
                    f"{a.xx:8.3f}{a.xy:8.3f}{a.xz:8.3f}"
                    f"  1.00  0.00          {a.element_name:2s}  \n")

        # Write CONECT for NZ-C bond
        if nz_atom_num and lnk_c_atom_num:
            f.write(f"CONECT{nz_atom_num:5d}{lnk_c_atom_num:5d}\n")

        f.write("TER\nEND\n")

    print(f"  {out_pdb}: {atom_num} atoms, CONECT {nz_atom_num}-{lnk_c_atom_num}")
    return out_pdb

def build_tleap(vname, combined_pdb):
    """Generate tleap input."""
    leap_in = f"{TLEAP}/build_{vname}.in"
    prmtop = f"{TLEAP}/{vname}.prmtop"
    inpcrd = f"{TLEAP}/{vname}.inpcrd"

    with open(leap_in, 'w') as f:
        f.write(f"""# tleap for exp-D: {vname}
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p

# Load combined PDB (with CONECT for NZ-C bond)
sys = loadPdb combined_{vname}.pdb

# Solvate
solvateOct sys TIP3PBOX 10.0
addIonsRand sys Na+ 0 Cl- 0

check sys

saveAmberParm sys {vname}.prmtop {vname}.inpcrd
savePdb sys {vname}_solvated.pdb
quit
""")
    return leap_in

if __name__ == "__main__":
    os.makedirs(TLEAP, exist_ok=True)
    for vname in VARIANTS:
        pdb = build_combined_pdb(vname)
        if pdb:
            # Copy PDB to tleap dir
            fname = os.path.basename(pdb)
            leap = build_tleap(vname, fname)
            print(f"  Tleap: {leap}")
    print("\nDone! Run tleap from TLEAP dir for each variant.")
