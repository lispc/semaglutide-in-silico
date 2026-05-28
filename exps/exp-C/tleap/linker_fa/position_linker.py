#!/usr/bin/env python3
"""Position linker-C18 diacid at HSA FA3 site and build tleap topology."""
import numpy as np, os
TL_DIR = "/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap"
OUT_DIR = f"{TL_DIR}/linker_fa"
os.chdir(OUT_DIR)

# Read mol2
with open(f"{TL_DIR}/linker_c18.mol2") as f:
    lines = f.readlines()

atom_start = None; atom_end = None
for i, line in enumerate(lines):
    if "@<TRIPOS>ATOM" in line: atom_start = i + 1
    if "@<TRIPOS>BOND" in line and atom_start: atom_end = i; break

atoms = []
for i in range(atom_start, atom_end):
    if not lines[i].strip(): continue
    p = lines[i].split()
    if len(p) >= 6:
        atoms.append({'idx': i, 'aid': int(p[0]), 'name': p[1],
                      'coord': np.array([float(p[2]), float(p[3]), float(p[4])]),
                      'type': p[5]})

# MYR 1003 FA3 template
myr_o1 = np.array([13.4, 8.3, 10.2]); myr_o2 = np.array([14.2, 10.0, 9.0])
myr_c1 = np.array([13.9, 9.5, 10.1]); myr_c14 = np.array([9.6, 12.6, 17.6])
carb_center = (myr_o1 + myr_o2 + myr_c1) / 3
chain_dir = (myr_c14 - myr_c1) / np.linalg.norm(myr_c14 - myr_c1)

# Find the distal carboxyl carbon and oxygens in linker-C18
# The distal COOH should be at one end of the molecule (furthest from the linker)
# Look for the carboxyl C (type 'c') that has two O neighbors (type 'o' or 'oh')
# The molecule has: α-COOH (on γGlu) and distal-COOH (on C18 chain)
# The distal one should be at the far end of the extended chain

# Simple approach: find all carboxyl carbons (type 'c'), pick the one furthest from molecule center
carb_c_atoms = [a for a in atoms if a['type'] == 'c']
mol_center = np.mean([a['coord'] for a in atoms], axis=0)

# The two carboxyl carbons: one near linker (γGlu α-COOH), one at distal end
# Pick the one furthest from center as the distal COOH
distal_c = max(carb_c_atoms, key=lambda a: np.linalg.norm(a['coord'] - mol_center))
# Find its O neighbors (within 1.5 Å)
distal_oxy = [a for a in atoms if a['type'] in ('o','oh','os') and
              np.linalg.norm(a['coord'] - distal_c['coord']) < 1.6]

print(f"Distal carboxyl C: {distal_c['name']} at {distal_c['coord']}")
print(f"Distal carboxyl O atoms: {len(distal_oxy)}")

if distal_oxy:
    fa_center = (distal_c['coord'] + sum(o['coord'] for o in distal_oxy)) / (1 + len(distal_oxy))
else:
    fa_center = distal_c['coord']

# Find approximate FA chain direction: from distal COOH towards the rest of the molecule
# Use the direction from distal C to the next C along the chain
other_coords = np.array([a['coord'] for a in atoms if a['aid'] != distal_c['aid']])
fa_dir = (np.mean(other_coords, axis=0) - distal_c['coord'])
fa_dir = fa_dir / np.linalg.norm(fa_dir)

# Rotation: align fa_dir to chain_dir
v = np.cross(fa_dir, chain_dir); s = np.linalg.norm(v); c = np.dot(fa_dir, chain_dir)
if s > 1e-8:
    vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    R = np.eye(3) + vx + vx@vx*((1-c)/(s*s))
else:
    R = np.eye(3)

# Apply transformation
translation = carb_center - fa_center
for a in atoms:
    a['coord'] = (a['coord'] - fa_center) @ R.T + carb_center
    # Update mol2 line
    parts = lines[a['idx']].split()
    lines[a['idx']] = f"{a['aid']:6d} {a['name']:6s} {a['coord'][0]:10.4f} {a['coord'][1]:10.4f} {a['coord'][2]:10.4f} {a['type']:4s} 1 LFA {float(parts[-1]):10.4f}\n"

# Save repositioned mol2
with open("linker_c18_fa3.mol2", 'w') as f:
    f.writelines(lines)

# Check min distance to HSA atoms
hsa_xyz = []
with open(f"{TL_DIR}/hsa_no_myr.pdb") as f:
    for line in f:
        if line.startswith("ATOM"):
            hsa_xyz.append(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))

lfa_xyz = [a['coord'] for a in atoms]
min_d = min(np.linalg.norm(l-h) for l in lfa_xyz[::5] for h in hsa_xyz[::10])
print(f"Min LFA-HSA distance: {min_d*10:.1f} Å")

# Build combined PDB
hsa_lines = []
with open(f"{TL_DIR}/hsa_no_myr.pdb") as f:
    for line in f:
        if line.startswith("ATOM"): hsa_lines.append(line)

with open("complex_linker.pdb", 'w') as f:
    f.write("REMARK HSA + linker-C18 diacid at FA3\n")
    aid = 0
    for line in hsa_lines:
        aid += 1; f.write(f"ATOM  {aid:5d}{line[11:]}")
    f.write("TER\n")
    for a in atoms:
        aid += 1; elem = a['name'][0] if a['name'][0] in 'CON' else 'C'
        f.write(f"HETATM{aid:5d} {a['name']:4s} LFA B   1    "
                f"{a['coord'][0]:8.3f}{a['coord'][1]:8.3f}{a['coord'][2]:8.3f}"
                f"  1.00  0.00          {elem:2s}  \n")
    f.write("TER\nEND\n")

print(f"Combined PDB: {aid} atoms total")
print("Saved: linker_c18_fa3.mol2, complex_linker.pdb")
