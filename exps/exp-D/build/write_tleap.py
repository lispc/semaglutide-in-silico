#!/usr/bin/env python3
"""Write tleap build scripts for each LYA variant."""
import os

TLEAP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
BUILD = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/build"

VARIANTS = ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

def write_tleap(vname):
    stripped_mol2 = f"{BUILD}/{vname}_stripped.mol2"
    if not os.path.exists(stripped_mol2):
        print(f"SKIP {vname}: {stripped_mol2} not found")
        return

    leap_in = f"{TLEAP}/build_{vname}.in"
    prmtop = f"{TLEAP}/{vname}.prmtop"
    inpcrd = f"{TLEAP}/{vname}.inpcrd"

    with open(leap_in, 'w') as f:
        f.write(f"""# tleap for exp-D: {vname}
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p

# Load ECD-peptide complex
complex = loadPdb ecd_pep.pdb

# Load stripped linker mol2
LNK = loadMol2 {vname}_stripped.mol2

# Combine into one system
sys = combine {{ complex LNK }}

# Remove HZ atoms from Lys26 (residue number 26)
# In combined unit, Lys26 is still residue 26 (or find by name)
# PDB numbering preserved: :26 refers to Lys26
remove sys sys:26@HZ1
remove sys sys:26@HZ2
remove sys sys:26@HZ3

# Bond Lys26 NZ to LNK carbonyl C (atom C in LNK residue)
# LNK is the last residue - use :LNK to reference it
bond sys:26@NZ sys:LNK@C

# Solvate
solvateOct sys TIP3PBOX 10.0
addIonsRand sys Na+ 0 Cl- 0

# Check
check sys

# Save
saveAmberParm sys {vname}.prmtop {vname}.inpcrd
savePdb sys {vname}.pdb
quit
""")
    print(f"Wrote: {leap_in}")

if __name__ == "__main__":
    os.makedirs(TLEAP, exist_ok=True)
    # Copy stripped mol2 files to tleap dir
    for vname in VARIANTS:
        write_tleap(vname)
        src = f"{BUILD}/{vname}_stripped.mol2"
        dst = f"{TLEAP}/{vname}_stripped.mol2"
        if os.path.exists(src):
            import shutil
            shutil.copy2(src, dst)
            print(f"  Copied {vname}_stripped.mol2 to tleap/")
    print("Done!")
