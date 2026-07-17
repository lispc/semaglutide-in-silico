#!/usr/bin/env python3
"""Frame-0 validation of freshly built c18_monoacid / c18_diacid topologies.

Checks (per exps/exp-C/exp-log.md hard-won lessons):
  1. FA atom count vs formula (mono 58 / diacid 60, all-H)
  2. FA connectivity: no isolated atoms, carboxyl C valence 3
  3. HSA completeness: 582 residues (1E7G 3-584), 582 CA
  4. Distal carboxyl O1D/O2D -> anchoring ARG guanidinium N (expect ~2.5-3.5 A;
     canonical ARG348/ARG485 = tleap residue 346/483)
  5. Net charge ~0; Na+/Cl- counts
Writes a report to frame0_validation.txt
"""
import numpy as np
import parmed as pmd

EXPECT = {'c18_monoacid': 58, 'c18_diacid': 60}
out = []

def rep(s=''):
    print(s)
    out.append(s)

for name, n_fa_expect in EXPECT.items():
    rep(f"{'='*64}\n{name}\n{'='*64}")
    prm = pmd.load_file(f"{name}.prmtop", f"{name}.inpcrd")
    rep(f"total atoms: {len(prm.atoms)}")

    # --- FA checks ---
    fa_res = next(r for r in prm.residues if r.name == 'FAH')
    fa_atoms = fa_res.atoms
    rep(f"FAH atoms: {len(fa_atoms)} (expect {n_fa_expect}) "
        f"-> {'OK' if len(fa_atoms) == n_fa_expect else 'FAIL'}")
    q_fa = sum(a.charge for a in fa_atoms)
    rep(f"FAH charge: {q_fa:+.3f}")

    isolated = [a.name for a in fa_atoms if len(a.bond_partners) == 0]
    rep(f"FAH isolated atoms: {isolated if isolated else 'none'}")
    for cname in ('C1X', 'C1D', 'C1P'):
        c = next((a for a in fa_atoms if a.name == cname), None)
        if c is not None:
            rep(f"  {cname} valence: {len(c.bond_partners)} "
                f"({sorted(p.name for p in c.bond_partners)}) -> "
                f"{'OK' if len(c.bond_partners) == 3 else 'FAIL'}")
    # every H bonded to a carbon
    bad_h = [a.name for a in fa_atoms if a.element == 1 and
             not any(p.element == 6 for p in a.bond_partners)]
    rep(f"FAH H not bonded to C: {bad_h if bad_h else 'none'}")
    # terminal connectivity: walk C1D/C1X -> far end
    c1 = next(a for a in fa_atoms if a.name in ('C1X', 'C1D'))
    seen, stack = set(), [c1]
    while stack:
        a = stack.pop()
        if a.idx in seen: continue
        seen.add(a.idx)
        stack.extend(p for p in a.bond_partners if p.residue.name == 'FAH')
    rep(f"FAH connected component size: {len(seen)} / {len(fa_atoms)} "
        f"-> {'OK' if len(seen) == len(fa_atoms) else 'FAIL'}")

    # --- HSA checks ---
    prot_res = [r for r in prm.residues if r.name not in
                ('FAH', 'WAT', 'Na+', 'Cl-', 'NA', 'CL', 'Na', 'Cl') and
                len(r.atoms) > 4]
    prot_res = [r for r in prot_res if any(a.name == 'CA' for a in r.atoms)]
    n_ca = sum(1 for r in prot_res for a in r.atoms if a.name == 'CA')
    rep(f"HSA residues: {len(prot_res)} (expect 582), CA count: {n_ca}")

    # --- ions / water ---
    from collections import Counter
    resn = Counter(r.name for r in prm.residues)
    rep(f"waters: {resn.get('WAT', 0)}, Na+: {resn.get('Na+', 0)}, Cl-: {resn.get('Cl-', 0)}")
    rep(f"net charge: {sum(a.charge for a in prm.atoms):+.4f}")
    box = prm.box
    rep(f"box: {box[0]:.1f} x {box[1]:.1f} x {box[2]:.1f} A")

    # --- carboxyl -> ARG distances ---
    xyz = prm.coordinates  # (n,3) Angstrom
    def coord(atom): return xyz[atom.idx]
    o_atoms = {a.name: a for a in fa_atoms if a.name in ('O1D', 'O2D')}
    arg_n = [(r, a) for r in prot_res if r.name == 'ARG'
             for a in r.atoms if a.name in ('NE', 'NH1', 'NH2')]
    best = None
    for oname, oa in sorted(o_atoms.items()):
        ds = sorted((np.linalg.norm(coord(ga) - coord(oa)), r.number + 1, ga.name)
                    for r, ga in arg_n)
        best = ds[0][0] if best is None else min(best, ds[0][0])
        top = ", ".join(f"ARG{rn}({nn}) {d:.2f}" for d, rn, nn in ds[:3])
        rep(f"  {oname} -> nearest ARG guanidinium N: {top}")
    ok = "OK" if best is not None and 2.0 <= best <= 3.5 else "FAIL"
    rep(f"  distal carboxyl anchor (best O..N): {best:.2f} A  [{ok}]")
    # tleap residue number = parmed index+1; canonical = tleap + 2 (PDB 3-584)
    rep("  (tleap ARG346 = canonical ARG348, tleap ARG483 = canonical ARG485"
        " = 'ARG482' in old 0-based-renumbered notes)")

    # --- min FA-HSA distance (clash check) ---
    fa_hvy = np.array([coord(a) for a in fa_atoms if a.element != 1])
    hsa_hvy = np.array([xyz[a.idx] for r in prot_res for a in r.atoms if a.element != 1])
    dmat = np.linalg.norm(fa_hvy[:, None, :] - hsa_hvy[None, :, :], axis=-1)
    rep(f"min FA-heavy..HSA-heavy distance: {dmat.min():.2f} A (clash if < 1.5)")
    rep()

with open('frame0_validation.txt', 'w') as f:
    f.write("\n".join(out) + "\n")
print("saved frame0_validation.txt")
