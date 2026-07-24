#!/usr/bin/env python3
"""Assemble the exp-G ternary complex solute PDB.

Applies the (NZ candidate, rotation) pair chosen by the joint docking
search in place_lnk_fa3.py (stored in tleap/ecd_dock_transform.npz) to the
rigid ECD-peptide template, merges with HSA, and writes the combined PDB
via ParmEd (so tleap's atom-mask parser reads it correctly).

Output: exps/exp-G/tleap/complex_ecd_pep_hsa.pdb + build_complex.in
tleap numbering (contiguous): HSA 1-582, ECD 583-682, peptide 683-708
(Lys26 = 699), LNK = 709.
"""
import os
import numpy as np
import parmed as pmd

HERE = os.path.dirname(os.path.abspath(__file__))
EXP_G = os.path.normpath(os.path.join(HERE, ".."))
TLEAP = os.path.join(EXP_G, "tleap")
HSA_PDB = "/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap/hsa_no_myr.pdb"
ECD_PEP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap/ecd_pep_nocap.pdb"

def rot_axis(axis, theta):
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * K @ K

dock = np.load(os.path.join(TLEAP, "ecd_dock_transform.npz"))
nz_target, axis, R0, deg = dock["nz_target"], dock["axis"], dock["R0"], float(dock["deg"])
print(f"dock transform: rot={deg:.0f} deg about axis {np.round(axis, 2)}")

hsa = pmd.load_file(HSA_PDB)
ecd = pmd.load_file(ECD_PEP)

nz0 = None
for a in ecd.atoms:
    if a.residue.chain == "B" and a.residue.number == 117 and a.name == "NZ":
        nz0 = np.array([a.xx, a.xy, a.xz]); break
assert nz0 is not None

e0 = np.array([[a.xx, a.xy, a.xz] for a in ecd.atoms])
e1 = (R0 @ (e0 - nz0).T).T + nz_target
R = rot_axis(axis, np.deg2rad(deg))
e2 = (R @ (e1 - nz_target).T).T + nz_target
for i, a in enumerate(ecd.atoms):
    a.xx, a.xy, a.xz = e2[i]

n_ecd = n_pep = 0
for r in ecd.residues:
    if r.chain == "A":
        r.chain = "R"; r.number = 601 + n_ecd; n_ecd += 1
    else:
        r.chain = "P"; r.number = 701 + n_pep; n_pep += 1
print(f"ECD residues: {n_ecd}, peptide residues: {n_pep}")

combined = hsa + ecd
out_pdb = os.path.join(TLEAP, "complex_ecd_pep_hsa.pdb")
combined.save(out_pdb, overwrite=True)
print(f"wrote {out_pdb}")

# sanity: Lys26 NZ position after transform
for a in combined.atoms:
    if a.residue.chain == "P" and a.residue.number == 717 and a.name == "NZ":
        d = np.linalg.norm(np.array([a.xx, a.xy, a.xz]) - nz_target)
        print(f"Lys26 NZ vs target: {d:.3f} A")
        assert d < 0.01

n_hsa_res = len(hsa.residues)
lys_tleap = n_hsa_res + n_ecd + 17
lnk_tleap = n_hsa_res + n_ecd + n_pep + 1
print(f"tleap numbering: Lys26 = {lys_tleap}, LNK = {lnk_tleap}")
assert lys_tleap == 699 and lnk_tleap == 709

with open(os.path.join(TLEAP, "build_complex.in"), "w") as f:
    f.write(f"""# tleap for exp-G ternary complex (pass 1: solvate only, count waters)
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadAmberParams frcmod.ionsjc_tip3p
loadAmberParams lya_link.frcmod
loadAmberParams lya_link_c8.frcmod

LNK = loadMol2 lnk_2oeg_fa3.mol2
prot = loadPdb complex_ecd_pep_hsa.pdb
complex = combine {{ prot LNK }}

# deprotonate Lys26 NZ, form NZ-C11 amide bond
remove complex complex.{lys_tleap}.HZ1
remove complex complex.{lys_tleap}.HZ2
remove complex complex.{lys_tleap}.HZ3
bond complex.{lys_tleap}.NZ complex.{lnk_tleap}.C11

solvateBox complex TIP3PBOX 12.0
saveAmberParm complex complex_nosalt.prmtop complex_nosalt.inpcrd
quit
""")
print("wrote build_complex.in (pass 1)")
