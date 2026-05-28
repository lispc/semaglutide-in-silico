#!/usr/bin/env python3
"""Build combined HSA + C18 diacid PDB for tleap."""
import os, numpy as np
os.chdir("/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap")

with open("hsa_no_myr.pdb") as f:
    hsa_lines = [l for l in f if l.startswith("ATOM")]

fah_atoms = []
with open("c18_diacid_fa3.mol2") as f:
    in_atom = False
    for line in f:
        if "@<TRIPOS>ATOM" in line: in_atom = True; continue
        if "@<TRIPOS>BOND" in line: break
        if in_atom and line.strip():
            p = line.split()
            if len(p) >= 6:
                fah_atoms.append({'name': p[1], 'x': float(p[2]), 'y': float(p[3]),
                                  'z': float(p[4]), 'elem': p[1][0]})

with open("complex_diacid.pdb", 'w') as f:
    f.write("REMARK HSA + C18 diacid at FA3\n")
    aid = 0
    for line in hsa_lines:
        aid += 1; f.write(f"ATOM  {aid:5d}{line[11:]}")
    f.write("TER\n")
    for a in fah_atoms:
        aid += 1; f.write(f"HETATM{aid:5d} {a['name']:4s} FAH B   1    "
                         f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                         f"  1.00  0.00          {a['elem']:2s}  \n")
    f.write("TER\nEND\n")

hs = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in hsa_lines])
fs = np.array([[a['x'], a['y'], a['z']] for a in fah_atoms])
md = min(np.linalg.norm(fc-hc) for fc in fs[::3] for hc in hs[::10])
print(f"Combined PDB: {aid} atoms, min FAH-HSA dist: {md*10:.1f} A")
print("Saved: complex_diacid.pdb")
