#!/usr/bin/env python3
"""Prepare ECD-peptide-Linker complexes for exp-D.
Strategy: clean 3IOL, rename Lys34→Arg, write PDB for tleap assembly.
tleap will add missing ARG atoms and handle LYA bonding.
"""
import os
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_D = f"{REPO}/exps/exp-D"
TLEAP = f"{EXP_D}/tleap"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def clean_and_mutate():
    """Extract ECD+peptide from 3IOL, rename Lys34→Arg, remove waters."""
    pdb = pmd.load_file(f"{REPO}/exps/exp-B/structures/3IOL.pdb")

    # Remove waters and non-protein
    keep_idx = []
    for i, a in enumerate(pdb.atoms):
        if a.residue.name in ('HOH', 'WAT'):
            continue
        if a.residue.chain in ('A', 'B'):
            keep_idx.append(i)
    pdb = pdb[keep_idx]

    # Mutate Lys34→Arg: rename residue, remove Lys-specific atoms
    for residue in pdb.residues:
        if residue.chain == 'B' and residue.number == 34 and residue.name == 'LYS':
            print(f"Mutating LYS{residue.number} → ARG (chain B)")
            residue.name = 'ARG'
            # Remove Lys-specific atoms that ARG doesn't have
            # ARG has: NE, CZ, NH1, NH2 instead of CE, NZ, HZ1-3
            # Common atoms: N, CA, C, O, CB, CG, CD
            for a in list(residue.atoms):
                if a.name in ('CE', 'NZ', 'HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3'):
                    # Keep CE and NZ positions but rename
                    if a.name == 'CE':
                        a.name = 'NE'
                        a.type = 'NE'
                    elif a.name == 'NZ':
                        a.name = 'CZ'
                        a.type = 'CZ'
                    elif a.name in ('HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3'):
                        # Remove hydrogens - tleap adds correct ones
                        continue  # will be removed below

            # Remove renamed H atoms
            to_remove = []
            for a in residue.atoms:
                if a.name in ('HE2', 'HE3', 'HZ1', 'HZ2', 'HZ3'):
                    to_remove.append(a)
            for a in to_remove:
                residue.delete_atom(a)

    # Save
    os.makedirs(TLEAP, exist_ok=True)
    out_pdb = f"{TLEAP}/ecd_pep.pdb"
    pdb.save(out_pdb, overwrite=True)
    print(f"Saved: {out_pdb} ({len(pdb.atoms)} atoms)")

    # Print info
    pep_res = [r for r in pdb.residues if r.chain == 'B']
    print(f"Peptide: {', '.join(f'{r.name}{r.number}' for r in pep_res)}")
    ecd_res = [r for r in pdb.residues if r.chain == 'A']
    print(f"ECD: {len(ecd_res)} residues ({ecd_res[0].name}{ecd_res[0].number}-{ecd_res[-1].name}{ecd_res[-1].number})")

    # Find Lys26 NZ index for later bonding
    for a in pdb.atoms:
        if a.residue.chain == 'B' and a.residue.number == 26 and a.name == 'NZ':
            print(f"Lys26 NZ: atom index {a.idx}, residue idx {a.residue.idx}")
            break

    return out_pdb

if __name__ == "__main__":
    os.makedirs(TLEAP, exist_ok=True)
    clean_and_mutate()

    # Generate tleap input for each variant
    for vname in VARIANTS:
        bcc_mol2 = f"{EXP_D}/build/{vname}_bcc.mol2"
        leap_in = f"{TLEAP}/build_{vname}.in"
        prmtop = f"{TLEAP}/{vname}.prmtop"
        inpcrd = f"{TLEAP}/{vname}.inpcrd"

        with open(leap_in, 'w') as f:
            f.write(f"""# tleap input for exp-D: {vname}
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p

# Load ECD + peptide complex
complex = loadPdb ecd_pep.pdb

# Load LYA linker
LYA = loadMol2 {vname}_bcc.mol2

# Bond LYA amide C to Lys26 NZ
# Find Lys26 in chain B
bond complex.116.NZ LYA.1.N
bond complex.116.NZ LYA.1.C

# Remove H atoms that would clash
# Lys26 HZ1/HZ2/HZ3 are removed automatically by bond
# LYA amide N loses one H

# Solvate
solvateOct complex TIP3PBOX 10.0
addIonsRand complex Na+ 0 Cl- 0

saveAmberParm complex {vname}.prmtop {vname}.inpcrd
quit
""")
        print(f"Wrote: {leap_in}")
        print(f"  NOTE: Bond command needs correct LYA atom index - verify after mol2 inspection")
