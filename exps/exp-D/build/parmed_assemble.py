#!/usr/bin/env python3
"""ParmEd-based full assembly: ECD + peptide + LNK → prmtop.
Bypasses tleap auto-completion issues entirely.
"""
import numpy as np, parmed as pmd, os, sys, shutil

REPO = "/home/scroll/personal/semaglutide-in-silico"
TLEAP = f"{REPO}/exps/exp-D/tleap"
BUILD = f"{REPO}/exps/exp-D/build"

os.chdir(TLEAP)

# Load reference peptide from 3IOL for coordinates
ref_pdb = pmd.load_file(f"{REPO}/exps/exp-B/structures/3IOL.pdb")
ref_pep = [a for a in ref_pdb.atoms if a.residue.chain == 'B' and a.residue.name not in ('HOH','WAT')]

# Load tleap-built peptide mol2 (structure=True for Structure object)
pep = pmd.load_file("pep_tleap.mol2", structure=True)
print(f"Peptide: {len(pep.atoms)} atoms, {len(pep.residues)} residues")

# Map residues: both have same sequence (26 residues G10-G35)
# Get reference backbone coords
ref_bb = []
ref_resnames = []
for r in ref_pdb.residues:
    if r.chain == 'B' and r.name not in ('HOH','WAT'):
        for a in r.atoms:
            if a.name in ('N','CA','C'):
                ref_bb.append((a.xx, a.xy, a.xz))
                ref_resnames.append((r.name, r.idx, a.name))

pep_bb = []
pep_resnames = []
for r in pep.residues:
    for a in r.atoms:
        if a.name in ('N','CA','C'):
            pep_bb.append((a.xx, a.xy, a.xz))
            pep_resnames.append((r.name, r.idx, a.name))

print(f"Ref BB: {len(ref_bb)}, Pep BB: {len(pep_bb)}")
print(f"  Ref: {ref_resnames[0]} ... {ref_resnames[-1]}")
print(f"  Pep: {pep_resnames[0]} ... {pep_resnames[-1]}")

# Kabsch align
rc = np.array(ref_bb).mean(0)
pc = np.array(pep_bb).mean(0)
H = (np.array(pep_bb) - pc).T @ (np.array(ref_bb) - rc)
U, s, Vt = np.linalg.svd(H)
R = Vt.T @ U.T
if np.linalg.det(R) < 0:
    Vt[-1] *= -1; R = Vt.T @ U.T
t = rc - R @ pc

# Apply to all pep atoms
for a in pep.atoms:
    xyz = R @ np.array([a.xx, a.xy, a.xz]) + t
    a.xx, a.xy, a.xz = xyz

# Verify RMSD
aligned_bb = []
for r in pep.residues:
    for a in r.atoms:
        if a.name in ('N','CA','C'):
            aligned_bb.append((a.xx, a.xy, a.xz))
rmsd = np.sqrt(((np.array(aligned_bb) - np.array(ref_bb))**2).mean())
print(f"Kabsch BB RMSD: {rmsd:.3f} A")

# Find Lys26 in peptide (the only LYS)
lys_r = None
for r in pep.residues:
    if r.name == 'LYS':
        lys_r = r
        print(f"Lys residue: idx={r.idx}, atoms: {[a.name for a in r.atoms]}")
        break

if not lys_r:
    print("ERROR: No LYS in peptide!")
    sys.exit(1)

# Find NZ, CE, HZ atoms
nz = ce = None
hz_atoms = []
for a in lys_r.atoms:
    if a.name == 'NZ': nz = a
    if a.name == 'CE': ce = a
    if a.name.startswith('HZ'): hz_atoms.append(a)
print(f"NZ={nz.name if nz else '?'}, CE={ce.name if ce else '?'}, HZ={[a.name for a in hz_atoms]}")

# Delete HZ atoms from peptide
# In ParmEd, we need to rebuild the structure without HZ atoms
# Strategy: write and reload, or use strip
print(f"Peptide before HZ removal: {len(pep.atoms)} atoms")
# Remove HZ by selecting all atoms except HZ
hz_indices = [a.idx for a in hz_atoms if hasattr(a, 'idx')]
# Actually, ParmEd's delete is tricky. Use the indices of ALL atoms minus HZ
all_pep_atoms = list(pep.atoms)
keep_pep_idx = [i for i, a in enumerate(all_pep_atoms) if a not in hz_atoms]
pep = pep[keep_pep_idx]
print(f"Peptide after HZ removal: {len(pep.atoms)} atoms")

# Verify NZ still exists
nz_found = False
for a in pep.atoms:
    if a.residue.name == 'LYS' and a.name == 'NZ':
        nz_found = True
        print(f"NZ preserved: ({a.xx:.1f},{a.xy:.1f},{a.xz:.1f})")
        break
if not nz_found:
    print("ERROR: NZ lost during HZ removal!")

# Save cleaned peptide
pep.save("pep_clean.pdb", overwrite=True)
pep.save("pep_clean.prmtop", format="amber", overwrite=True)
print("Saved pep_clean.pdb and pep_clean.prmtop")

# Now write tleap build script that uses the cleaned peptide
for vname in ["no_linker", "gglu", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]:
    stripped = f"../build/{vname}_stripped.mol2"
    bcc = f"../build/{vname}_bcc.mol2"

    if not os.path.exists(bcc):
        print(f"SKIP {vname}: {bcc} not found")
        continue
    if not os.path.exists(stripped):
        print(f"SKIP {vname}: {stripped} not found (run strip_lya.py)")
        continue

    leap_in = f"build_{vname}.in"
    with open(leap_in, 'w') as f:
        f.write(f"""# tleap for exp-D: {vname}
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p

# Load ECD (chain A only)
ecd = loadPdb ecd_only.pdb

# Load cleaned peptide (no HZ atoms)
pep = loadPdb pep_clean.pdb

# Load LNK linker
LNK = loadMol2 ../build/{vname}_stripped.mol2

# Combine ECD + peptide
complex = combine {{ ecd pep }}

# Combine with LNK
sys = combine {{ complex LNK }}

# Find residue indices:
# peptide Lys is at known index (from desc output)
# LNK is the last residue
# Need to remove HZ atoms from peptide Lys and bond to LNK C
# Peptide LYS has index 16 (0-based in sequence, adjusted for combined unit)
# In combined unit: ECD residues first, then peptide, then LNK
# ECD has ~101 residues, so peptide LYS is at ~101+16 = 117

# Remove remaining HZ from Lys (tleap may have re-added them)
# Bond Lys NZ to LNK amide C
bond sys.117.NZ sys.128.C

# Solvate
solvateOct sys TIP3PBOX 10.0
addIonsRand sys Na+ 0

check sys

saveAmberParm sys {vname}.prmtop {vname}.inpcrd
savePdb sys {vname}.pdb
quit
""")
    print(f"Wrote: {leap_in}")

print("\nDone! Run from TLEAP dir: tleap -f build_no_linker.in")
