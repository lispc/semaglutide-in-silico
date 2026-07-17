#!/usr/bin/env python3
"""Prepare hsa_no_myr.pdb from 1E7G for tleap.

- Keep protein ATOM records, chain A (residues 3-584, 582 residues)
- Drop all HETATM (8x MYR, HOH) -- C18 FA is rebuilt separately
- Rename disulfide-bonded CYS -> CYX (from SSBOND records) so tleap
  does not expect HG on bonded cysteines (ff14SB)
- No H atoms in 1E7G; tleap adds them from templates
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "structures", "1E7G.pdb")
DST = os.path.join(HERE, "hsa_no_myr.pdb")

# Parse SSBOND records -> set of cysteine resids in disulfides
ss_cys = set()
with open(SRC) as f:
    for line in f:
        if line.startswith("SSBOND"):
            p = line.split()  # SSBOND n CYS A r1 CYS A r2 ...
            ss_cys.add(int(p[4]))
            ss_cys.add(int(p[7]))

out = []
n_cyx = 0
with open(SRC) as f:
    for line in f:
        if not line.startswith("ATOM"):
            continue
        resn = line[17:20]
        resi = int(line[22:26])
        if resn == "CYS" and resi in ss_cys:
            line = line[:17] + "CYX" + line[20:]
            n_cyx += 1
        out.append(line)
out.append("TER\n")
out.append("END\n")

with open(DST, "w") as f:
    f.writelines(out)

resids = sorted({int(l[22:26]) for l in out if l.startswith("ATOM")})
print(f"Wrote {DST}")
print(f"  ATOM lines: {sum(1 for l in out if l.startswith('ATOM'))}")
print(f"  residues: {resids[0]}..{resids[-1]} ({len(resids)} total)")
print(f"  CYX atoms renamed: {n_cyx} (disulfide CYS: {sorted(ss_cys)})")
