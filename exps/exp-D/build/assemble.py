#!/usr/bin/env python3
"""Assemble ECD-peptide + LYA linker complexes using ParmEd, solvate with tleap.
Bonds peptide Lys26 CE → LYA N10 (amide N), removes HZ atoms.
"""
import os, sys, subprocess, shutil
import parmed as pmd
import numpy as np

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_D = f"{REPO}/exps/exp-D"
TLEAP = f"{EXP_D}/tleap"
BUILD = f"{EXP_D}/build"
MD = f"{EXP_D}/md"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def assemble(variant_name):
    """Assemble ECD-pep + LYA into a combined PDB for tleap solvation."""
    lya_mol2 = f"{BUILD}/{variant_name}_bcc.mol2"
    if not os.path.exists(lya_mol2):
        print(f"  SKIP: {lya_mol2} not found")
        return None

    print(f"\n=== {variant_name} ===")

    # Load ECD-peptide
    complex_pdb = f"{TLEAP}/ecd_pep.pdb"
    comp = pmd.load_file(complex_pdb)

    # Load LYA
    lya = pmd.load_file(lya_mol2)

    # Find Lys26 atoms in complex
    lys26 = None
    for r in comp.residues:
        if r.chain == 'B' and r.number == 26 and r.name == 'LYS':
            lys26 = r
            break
    if not lys26:
        print("  ERROR: Lys26 not found")
        return None

    # Find key atoms
    nz = ce = None
    hz_atoms = []
    for a in lys26.atoms:
        if a.name == 'NZ': nz = a
        if a.name == 'CE': ce = a
        if a.name in ('HZ1', 'HZ2', 'HZ3'): hz_atoms.append(a)
    print(f"  Lys26 NZ={nz.idx if nz else '?'}, CE={ce.idx if ce else '?'}, HZ={len(hz_atoms)}")

    # Find LYA atoms: N10 = NZ-equivalent, C11 = amide C
    # Atoms C1-C9 = ACE+backbone+sidechain, keep N10 onward (NZ+linker+C18+NME)
    lya_n10 = lya_c11 = lya_nme_n = lya_nme_c = None
    for a in lya.atoms:
        if a.name == 'N10': lya_n10 = a
        if a.name == 'C11': lya_c11 = a
        if a.name == 'N34': lya_nme_n = a
        if a.name == 'C35': lya_nme_c = a
    print(f"  LYA N10={lya_n10.idx if lya_n10 else '?'}, C11={lya_c11.idx if lya_c11 else '?'}")

    if not all([nz, ce, lya_n10, lya_c11]):
        print("  ERROR: missing key atoms")
        return None

    # Strategy: remove HZ from Lys26, bond CE to LYA N10
    # LYA N10 already has the -C11(=O)-linker-C18 tail attached
    # Final connectivity: ...Lys26 CE - N10 - C11(=O) - linker - C18 - NME

    # Keep atoms: complex (minus HZ) + LYA (N10 through end, minus ACE/backbone/sidechain C1-C9)
    # Also remove one H from LYA N10 (it has 1-2 H atoms)

    # Which H atoms on N10 to remove?
    h_on_n10 = [a for a in lya.atoms if a.element_name == 'H' and any(b is lya_n10 for b in a.bond_partners)]
    print(f"  H on N10: {[a.name for a in h_on_n10]}")

    # Complex atoms: exclude HZ
    comp_keep_idx = [i for i in range(len(comp.atoms))
                     if comp.atoms[i] not in hz_atoms]

    # LYA atoms to keep: N10 onward (idx >= lya_n10.idx), minus one H on N10
    # (The H on N10 is replaced by the bond to CE)
    lya_start = lya_n10.idx
    lya_keep_idx = []
    for i in range(lya_start, len(lya.atoms)):
        a = lya.atoms[i]
        # Skip one H on N10 (keep the other if there are 2)
        if a in h_on_n10[:1]:  # skip first H
            continue
        lya_keep_idx.append(i)

    # Build combined atom list
    combined_atoms = []
    # ECD + peptide (minus HZ)
    for i in comp_keep_idx:
        a = comp.atoms[i]
        combined_atoms.append({
            'name': a.name, 'elem': a.element_name,
            'resname': a.residue.name, 'resid': a.residue.idx + 1,
            'chain': a.residue.chain,
            'x': a.xx, 'y': a.xy, 'z': a.xz,
            'charge': a.charge if hasattr(a, 'charge') else 0.0,
        })

    # LYA atoms (N10 onward, minus one H)
    # Rename N10 to NZ for consistency with peptide naming
    lya_offset = len(combined_atoms)
    for i in lya_keep_idx:
        a = lya.atoms[i]
        name = a.name
        resname = 'LYA'
        if name == 'N10':
            name = 'NZ'  # rename to peptide naming convention
        combined_atoms.append({
            'name': name, 'elem': a.element_name,
            'resname': resname, 'resid': 27,  # after peptide, before ECD
            'chain': 'B',
            'x': a.xx, 'y': a.xy, 'z': a.xz,
            'charge': a.charge if hasattr(a, 'charge') else 0.0,
        })

    print(f"  Combined: {len(combined_atoms)} atoms ({len(comp_keep_idx)} complex + {len(lya_keep_idx)} LYA)")

    # Write PDB
    out_pdb = f"{TLEAP}/combined_{variant_name}.pdb"
    with open(out_pdb, 'w') as f:
        f.write(f"REMARK ECD-peptide + {variant_name} linker\n")
        for i, a in enumerate(combined_atoms):
            f.write(f"ATOM  {i+1:5d} {a['name']:4s} {a['resname']:3s} "
                    f"{a['chain']}{a['resid']:4d}    "
                    f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                    f"  1.00  0.00          {a['elem']:2s}  \n")
        # Write CONECT record for CE-NZ bond
        ce_atom_num = None
        nz_atom_num = None
        for i, a in enumerate(combined_atoms):
            if a['name'] == 'CE' and a['resname'] == 'LYS':
                ce_atom_num = i + 1
            if a['name'] == 'NZ' and a['resname'] == 'LYA':
                nz_atom_num = i + 1
        if ce_atom_num and nz_atom_num:
            f.write(f"CONECT{ce_atom_num:5d}{nz_atom_num:5d}\n")
        f.write("TER\nEND\n")
    print(f"  Saved: {out_pdb}")
    print(f"  CE atom #{ce_atom_num}, NZ atom #{nz_atom_num}")

    return out_pdb

def build_tleap_script(variant_name, combined_pdb):
    """Generate tleap build script."""
    leap_in = f"{TLEAP}/build_{variant_name}.in"
    prmtop = f"{TLEAP}/{variant_name}.prmtop"
    inpcrd = f"{TLEAP}/{variant_name}.inpcrd"

    with open(leap_in, 'w') as f:
        f.write(f"""# tleap input for exp-D: {variant_name}
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p

# Load combined ECD-peptide-LYA complex
complex = loadPdb {combined_pdb}

# Load LYA frcmod (AM1-BCC charges)
# complex already has LYA atoms with correct charges
# Need to set LYA bond parameters from GAFF2
loadAmberParams {BUILD}/{variant_name}_bcc.frcmod

# Solvate
solvateOct complex TIP3PBOX 10.0
addIonsRand complex Na+ 0 Cl- 0

# Check
check complex

saveAmberParm complex {variant_name}.prmtop {variant_name}.inpcrd
quit
""")
    print(f"  Wrote: {leap_in}")
    return leap_in

if __name__ == "__main__":
    os.makedirs(TLEAP, exist_ok=True)
    os.makedirs(MD, exist_ok=True)

    for vname in VARIANTS:
        pdb = assemble(vname)
        if pdb:
            # Copy LYA mol2 to tleap dir
            src = f"{BUILD}/{vname}_bcc.mol2"
            dst = f"{TLEAP}/{vname}_bcc.mol2"
            if os.path.exists(src):
                shutil.copy2(src, dst)

            # Generate frcmod (use default GAFF2 params)
            # tleap can load mol2 directly with loadMol2
            # build_tleap_script(vname, os.path.basename(pdb))

    print("\nDone! Next: run tleap for each variant")
