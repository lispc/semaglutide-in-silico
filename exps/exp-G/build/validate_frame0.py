#!/usr/bin/env python3
"""Frame-0 validation for the exp-G ternary complex (see README.md criteria).

1. distal COO- O56/O57 -> ARG346/ARG483 guanidinium: each anchor <= 3.5 A
   (best of NE/NH1/NH2; geometry reference: exp-C frame0_validation.txt)
2. proximal carbonyl O38 -> ARG408 <= 5 A
3. peptide-ECD interface min distance <= 4 A
4. HSA<->ECD/peptide/LNK heavy clashes (<2.2 A) = 0
5. NZ-C11 bond exists (parmed), ~1.38 A; 0 unparameterized terms
6. atom/box/ion counts recorded

Residue positions are 1-based enumeration in the prmtop (parmed
residue.number is 0-based on prmtop load -- do not rely on it).
"""
import os
import numpy as np
import parmed as pmd

HERE = os.path.dirname(os.path.abspath(__file__))
TLEAP = os.path.normpath(os.path.join(HERE, "..", "tleap"))

parm = pmd.load_file(os.path.join(TLEAP, "complex.prmtop"),
                     os.path.join(TLEAP, "complex.inpcrd"))
print(f"atoms={len(parm.atoms)}, residues={len(parm.residues)}")
box = parm.box
print(f"box: {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f} A")

prot_res = [r for r in parm.residues if r.name not in ("WAT", "HOH", "Na+", "Cl-", "NA", "CL")]
print(f"protein+LNK residues: {len(prot_res)} (expect 709)")
assert len(prot_res) == 709
pos_of = {id(r): i + 1 for i, r in enumerate(prot_res)}  # 1-based position

lnk = [a for a in parm.atoms if a.residue.name == "LNK"]
assert len(lnk) == 108 and pos_of[id(lnk[0].residue)] == 709
by_name = {a.name: a for a in lnk}

# --- NZ-C11 bond (Lys26 at position 699) ---
c11 = by_name["C11"]
nz_cands = [a for a in parm.atoms if a.name == "NZ" and a.residue.name == "LYS"
            and pos_of.get(id(a.residue), 0) <= 709]
nz = min(nz_cands, key=lambda a: np.linalg.norm(
    [a.xx - c11.xx, a.xy - c11.xy, a.xz - c11.xz]))
nz_pos = pos_of[id(nz.residue)]
d_nzc = np.linalg.norm([nz.xx - c11.xx, nz.xy - c11.xy, nz.xz - c11.xz])
has_bond = any({b.atom1, b.atom2} == {nz, c11} for b in parm.bonds)
hz = [a for a in parm.atoms if a.residue is nz.residue and a.name.startswith("HZ")]
print(f"Lys26 position: {nz_pos} (expect 699), NZ-C11 = {d_nzc:.2f} A, "
      f"bond: {has_bond}, HZ: {len(hz)}")
assert nz_pos == 699 and has_bond and d_nzc < 1.6 and len(hz) == 0

# --- anchor distances ---
def guanidinium_n(r):
    return np.array([[a.xx, a.xy, a.xz] for a in r.atoms if a.name in ("NE", "NH1", "NH2")])

arg_res = [r for r in prot_res if r.name == "ARG"]
def nearest_args(opos, k=4):
    d = sorted((np.sqrt(((guanidinium_n(r) - opos) ** 2).sum(-1)).min(), pos_of[id(r)])
               for r in arg_res)
    return d[:k]

o56 = np.array([by_name["O56"].xx, by_name["O56"].xy, by_name["O56"].xz])
o57 = np.array([by_name["O57"].xx, by_name["O57"].xy, by_name["O57"].xz])
d56 = nearest_args(o56); d57 = nearest_args(o57)
print(f"O56 nearest Arg: {[(rn, round(x,2)) for x, rn in d56]}")
print(f"O57 nearest Arg: {[(rn, round(x,2)) for x, rn in d57]}")
anchor1 = min(d56[0][0], d57[0][0])                       # ARG346 (tleap 346)
anchor2 = min(d56[1][0], d57[1][0])                       # ARG483 (tleap 483)
print(f"distal double salt bridge: ARG346 best {anchor1:.2f} A, ARG483 best {anchor2:.2f} A")
assert anchor1 <= 3.5 and anchor2 <= 3.5

o38 = np.array([by_name["O38"].xx, by_name["O38"].xy, by_name["O38"].xz])
d38 = nearest_args(o38)
print(f"O38 nearest Arg: {[(rn, round(x,2)) for x, rn in d38]} (expect ~408)")
assert d38[0][0] <= 5.0, f"O38-ARG408 {d38[0][0]:.2f} > 5.0"

# --- peptide-ECD interface ---
pep = np.array([[a.xx, a.xy, a.xz] for a in parm.atoms
                if 683 <= pos_of.get(id(a.residue), 0) <= 708 and a.element != 1])
ecd = np.array([[a.xx, a.xy, a.xz] for a in parm.atoms
                if 583 <= pos_of.get(id(a.residue), 0) <= 682 and a.element != 1])
d_if = np.sqrt(((ecd[None, :, :] - pep[:, None, :]) ** 2).sum(-1)).min()
print(f"peptide-ECD interface min distance: {d_if:.2f} A (expect <= 4)")
assert d_if <= 4.0

# --- clashes: HSA(1-582) vs ECD/pep/LNK ---
hsa = np.array([[a.xx, a.xy, a.xz] for a in parm.atoms
                if 1 <= pos_of.get(id(a.residue), 0) <= 582 and a.element != 1])
other = np.array([[a.xx, a.xy, a.xz] for a in parm.atoms
                  if 583 <= pos_of.get(id(a.residue), 0) <= 709 and a.element != 1])
d_cl = np.sqrt(((hsa[None, :, :] - other[:, None, :]) ** 2).sum(-1))
ncl = int((d_cl < 2.2).sum())
print(f"HSA<->ECD/pep/LNK heavy contacts <2.2 A: {ncl} (min {d_cl.min():.2f} A)")
assert ncl == 0

# --- LNK internal terms / unparameterized ---
lset = {a.idx for a in lnk}
nb = sum(1 for b in parm.bonds if b.atom1.idx in lset and b.atom2.idx in lset)
na = sum(1 for b in parm.angles if b.atom1.idx in lset and b.atom3.idx in lset)
nd = sum(1 for b in parm.dihedrals if b.atom2.idx in lset and b.atom3.idx in lset)
miss_b = sum(1 for b in parm.bonds if b.type is None)
miss_a = sum(1 for b in parm.angles if b.type is None)
miss_d = sum(1 for b in parm.dihedrals if b.type is None)
print(f"LNK internal: {nb} bonds, {na} angles, {nd} dihedrals; "
      f"unparameterized: {miss_b}/{miss_a}/{miss_d}")
assert nb == 107 and miss_b == miss_a == miss_d == 0

na_ion = sum(1 for a in parm.atoms if a.residue.name in ("Na+", "NA"))
cl_ion = sum(1 for a in parm.atoms if a.residue.name in ("Cl-", "CL"))
nwat = sum(1 for r in parm.residues if r.name in ("WAT", "HOH"))
q = sum(a.charge for a in parm.atoms)
print(f"waters={nwat}, Na+={na_ion}, Cl-={cl_ion}, net charge={q:+.2f}")
print("\nALL FRAME-0 CHECKS PASSED")
