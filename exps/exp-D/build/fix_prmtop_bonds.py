#!/usr/bin/env python3
"""Verify (and if needed, repair) exp-D prmtop files after tleap (2026-07-17).

Checks per variant (frame-0 validation):
1. NZ-C11 amide bond exists between Lys26 NZ and LNK carbonyl C (~1.4 A)
2. Lys26 has no HZ atoms (tleap `remove` has a documented silent-failure bug)
3. LNK internal bond/angle/dihedral counts are consistent
4. Total atom count ~36k; ion counts; net charge
5. The C8-N3-c-o / C8-N3-c-c3 dihedrals exist (tleap cannot parameterize
   them; patched here with the lya_link.frcmod values)
6. No steric clash among solute heavy atoms (excluding bonded/angle pairs)

If the tleap `bond`/`remove` silently failed, this script repairs:
adds NZ-C11 bond (BondType 427.0, 1.38 from lya_link.frcmod) and deletes HZ.
"""
import os, sys
import numpy as np
import parmed as pmd
from parmed.topologyobjects import Bond, BondType, DihedralType

TLEAP = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tleap"))
VARIANTS = ["no_linker", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]
EXPECTED_LNK = {"no_linker": 53, "gglu_1oeg": 87, "gglu_2oeg": 108, "gglu_3oeg": 129}

def verify(v):
    prmtop = os.path.join(TLEAP, f"{v}.prmtop")
    inpcrd = os.path.join(TLEAP, f"{v}.inpcrd")
    parm = pmd.load_file(prmtop, inpcrd)
    print(f"\n=== {v}: {len(parm.atoms)} atoms, {len(parm.residues)} residues ===")

    lnk = [a for a in parm.atoms if a.residue.name == "LNK"]
    assert len(lnk) == EXPECTED_LNK[v], f"LNK atom count {len(lnk)} != {EXPECTED_LNK[v]}"
    c11 = next(a for a in lnk if a.name == "C11")

    # NZ closest to C11 among all LYS
    cands = [a for a in parm.atoms if a.name == "NZ" and a.residue.name == "LYS"]
    nz = min(cands, key=lambda a: np.linalg.norm(
        np.array([a.xx - c11.xx, a.xy - c11.xy, a.xz - c11.xz])))
    d_nzc = np.linalg.norm(np.array([nz.xx - c11.xx, nz.xy - c11.xy, nz.xz - c11.xz]))
    print(f"  Lys26 NZ = residue {nz.residue.name} {nz.residue.number}, NZ-C11 = {d_nzc:.2f} A")

    # 1) NZ-C11 bond
    has_bond = any({b.atom1, b.atom2} == {nz, c11} for b in parm.bonds)
    if not has_bond:
        print("  REPAIR: tleap bond silently failed -> adding NZ-C11 BondType(427.0, 1.38)")
        parm.bonds.append(Bond(nz, c11, type=BondType(427.0, 1.38)))

    # 2) HZ on the acylated Lys
    hz = [a for a in parm.atoms if a.residue is nz.residue and a.name.startswith("HZ")]
    if hz:
        print(f"  REPAIR: tleap remove silently failed -> deleting {len(hz)} HZ atoms")
        drop = {a.idx for a in hz}
        parm = parm[[i for i in range(len(parm.atoms)) if i not in drop]]
        lnk = [a for a in parm.atoms if a.residue.name == "LNK"]
        c11 = next(a for a in lnk if a.name == "C11")
        nz = min([a for a in parm.atoms if a.name == "NZ" and a.residue.name == "LYS"],
                 key=lambda a: np.linalg.norm(np.array([a.xx - c11.xx, a.xy - c11.xy, a.xz - c11.xz])))
        has_bond = any({b.atom1, b.atom2} == {nz, c11} for b in parm.bonds)
        if not has_bond:
            parm.bonds.append(Bond(nz, c11, type=BondType(427.0, 1.38)))
    else:
        print("  HZ on Lys26: 0 (tleap remove OK)")

    # other lysines keep their HZ?
    n_other_lys_hz = sum(1 for a in parm.atoms
                         if a.residue.name == "LYS" and a.residue is not nz.residue
                         and a.name.startswith("HZ"))
    print(f"  HZ on other Lys residues: {n_other_lys_hz} (expect 6)")

    # 3) LNK internal bonds / angles / dihedrals
    lset = set(a.idx for a in lnk)
    nb = sum(1 for b in parm.bonds if b.atom1.idx in lset and b.atom2.idx in lset)
    na = sum(1 for b in parm.angles if b.atom1.idx in lset and b.atom3.idx in lset)
    nd = sum(1 for b in parm.dihedrals if b.atom2.idx in lset and b.atom3.idx in lset)
    print(f"  LNK internal: {nb} bonds, {na} angles, {nd} dihedrals (2-3 in LNK)")

    # 4) ions / charge
    na_ion = sum(1 for a in parm.atoms if a.residue.name in ("Na+", "NA"))
    cl_ion = sum(1 for a in parm.atoms if a.residue.name in ("Cl-", "CL"))
    nwat = sum(1 for r in parm.residues if r.name in ("WAT", "HOH"))
    q = sum(a.charge for a in parm.atoms)
    print(f"  waters={nwat}, Na+={na_ion}, Cl-={cl_ion}, net charge={q:+.3f}")

    # 5) C8-N3-c dihedrals: find dihedrals spanning NZ-C11, check/patch params
    patched = 0
    for d in parm.dihedrals:
        atoms = (d.atom1, d.atom2, d.atom3, d.atom4)
        names = (atoms[0].type, atoms[1].type, atoms[2].type, atoms[3].type)
        mid = {atoms[1], atoms[2]}
        if mid == {nz, c11} or (atoms[1] is nz and atoms[2] is c11) or (atoms[1] is c11 and atoms[2] is nz):
            if d.type is None:
                d.type = DihedralType(2.5, 2, 180.0)
                patched += 1
                print(f"  PATCH: dihedral {'-'.join(names)} had no params -> DihedralType(2.5, 2, 180)")
    # any other dihedral with no params?
    missing = [d for d in parm.dihedrals if d.type is None]
    if missing:
        print(f"  WARNING: {len(missing)} dihedrals without parameters:")
        for d in missing[:8]:
            print(f"    {d.atom1.type}-{d.atom2.type}-{d.atom3.type}-{d.atom4.type} "
                  f"({d.atom1.name}-{d.atom2.name}-{d.atom3.name}-{d.atom4.name})")
    # angles without params?
    ma = [a for a in parm.angles if a.type is None]
    mb = [b for b in parm.bonds if b.type is None]
    print(f"  unparameterized: {len(mb)} bonds, {len(ma)} angles, {len(missing)} dihedrals")

    # save if anything changed
    if not has_bond or hz or patched:
        parm.save(prmtop, format="amber", overwrite=True)
        parm.save(inpcrd, format="rst7", overwrite=True)
        print(f"  saved repaired {prmtop}")
        parm = pmd.load_file(prmtop, inpcrd)  # reload for clash check

    # 6) clash check: solute heavy atoms, exclude 1-2/1-3 pairs
    sol = [a for a in parm.atoms if a.residue.name not in ("WAT", "HOH", "Na+", "Cl-", "NA", "CL")
           and a.element != 1]
    pos = np.array([[a.xx, a.xy, a.xz] for a in sol])
    bonded = set()
    for b in parm.bonds:
        bonded.add(frozenset((b.atom1.idx, b.atom2.idx)))
    ang13 = set()
    for b in parm.angles:
        ang13.add(frozenset((b.atom1.idx, b.atom3.idx)))
    idx = {a.idx: k for k, a in enumerate(sol)}
    worst = []
    for i in range(len(sol)):
        d = np.linalg.norm(pos - pos[i], axis=1)
        for j in np.where((d < 1.9) & (d > 0))[0]:
            if j <= i: continue
            p = frozenset((sol[i].idx, sol[j].idx))
            if p in bonded or p in ang13: continue
            worst.append((d[j], sol[i], sol[j]))
    worst.sort()
    if worst:
        print(f"  close contacts (<1.9 A, non-bonded): {len(worst)}, worst:")
        for d, a1, a2 in worst[:5]:
            print(f"    {d:.2f} A  {a1.residue.name}{a1.residue.number}:{a1.name} - "
                  f"{a2.residue.name}{a2.residue.number}:{a2.name}")
    else:
        print("  no non-bonded heavy-atom contacts < 1.9 A")
    return True

if __name__ == "__main__":
    for v in VARIANTS:
        verify(v)
    print("\nDone.")
