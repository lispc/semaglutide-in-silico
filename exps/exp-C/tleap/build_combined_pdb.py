#!/usr/bin/env python3
"""Build combined HSA+FAH PDB in single coordinate frame."""
import numpy as np, os

os.chdir("/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap")

# Read HSA PDB
with open("hsa_no_myr.pdb") as f:
    hsa_lines = [l for l in f if l.startswith("ATOM")]

# Read FAH from positioned mol2
fah_atoms = []
with open("c18_monoacid_fa3.mol2") as f:
    in_atom = False
    for line in f:
        if "@<TRIPOS>ATOM" in line: in_atom = True; continue
        if "@<TRIPOS>BOND" in line: break
        if in_atom and line.strip():
            p = line.split()
            if len(p) >= 6:
                fah_atoms.append({'name': p[1], 'x': float(p[2]), 'y': float(p[3]),
                                  'z': float(p[4]), 'elem': p[1][0]})

# Write combined PDB
with open("complex_mono.pdb", 'w') as f:
    f.write("REMARK HSA + C18 monoacid at FA3\n")
    aid = 0
    for line in hsa_lines:
        aid += 1; f.write(f"ATOM  {aid:5d}{line[11:]}")
    f.write("TER\n")
    for a in fah_atoms:
        aid += 1; f.write(f"HETATM{aid:5d} {a['name']:4s} FAH B   1    "
                         f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                         f"  1.00  0.00          {a['elem']:2s}  \n")
    f.write("TER\nEND\n")

# Verify positions
import numpy as np
hsa_xyz = []; fah_xyz = []
for l in hsa_lines:
    hsa_xyz.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
for a in fah_atoms:
    fah_xyz.append(np.array([a['x'], a['y'], a['z']]))
hc = np.mean(hsa_xyz, axis=0); fc = np.mean(fah_xyz, axis=0)
min_d = min(np.linalg.norm(fc2 - hc2) for fc2 in fah_xyz[::3] for hc2 in hsa_xyz[::10])
print(f"HSA center: {hc}, FAH center: {fc}")
print(f"Distance: {np.linalg.norm(fc-hc):.1f} A, Min FAH-HSA: {min_d*10:.1f} A")
print("Combined PDB saved: complex_mono.pdb")
