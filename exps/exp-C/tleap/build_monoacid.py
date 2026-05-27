#!/usr/bin/env python3
"""Build C18 monoacid (stearate) mol2 at HSA FA3 site."""
import numpy as np

# MYR 1003 template at FA3
myr_o1 = np.array([13.4, 8.3, 10.2]); myr_o2 = np.array([14.2, 10.0, 9.0])
myr_c1 = np.array([13.9, 9.5, 10.1]); myr_c14 = np.array([9.6, 12.6, 17.6])
carboxyl_center = (myr_o1 + myr_o2 + myr_c1) / 3
chain_dir_u = (myr_c14 - myr_c1) / np.linalg.norm(myr_c14 - myr_c1)

n_c = 18; bond_len_cc = 1.54; angle_cc = 112.0 * np.pi / 180

# Build C chain at origin
c_pos = [np.array([0.0, 0.0, 0.0]), np.array([bond_len_cc, 0.0, 0.0])]
for i in range(2, n_c):
    prev_u = (c_pos[-1] - c_pos[-2]) / np.linalg.norm(c_pos[-1] - c_pos[-2])
    rot_axis = np.array([0.0, 0.0, 1.0])
    new_dir = np.cos(angle_cc) * prev_u + np.sin(angle_cc) * np.cross(rot_axis, prev_u)
    c_pos.append(c_pos[-1] + new_dir * bond_len_cc)

# Carboxyl at C1 end (distal)
c1 = c_pos[0]; c2 = c_pos[1]
to_c2_u = (c2 - c1) / np.linalg.norm(c2 - c1)
perp = np.cross(to_c2_u, np.array([0.,0.,1.]))
if np.linalg.norm(perp) < 0.01: perp = np.cross(to_c2_u, np.array([1.,0.,0.]))
perp = perp / np.linalg.norm(perp)

ccx_pos = c1 - to_c2_u * 1.52
o1d_pos = ccx_pos + perp * 1.08
o2d_pos = ccx_pos - perp * 1.08

print(f"O1D-O2D distance: {np.linalg.norm(o1d_pos-o2d_pos):.2f} A")

# Build atom list (heavy atoms only — tleap adds H from GAFF2)
atoms = []
aid = 0
atoms.append({'aid': (aid:=aid+1), 'name': 'C1X', 'coord': ccx_pos, 'type': 'c',  'charge': 0.70})
atoms.append({'aid': (aid:=aid+1), 'name': 'O1D', 'coord': o1d_pos, 'type': 'o',  'charge': -0.80})
atoms.append({'aid': (aid:=aid+1), 'name': 'O2D', 'coord': o2d_pos, 'type': 'o',  'charge': -0.80})
for i, cp in enumerate(c_pos):
    atoms.append({'aid': (aid:=aid+1), 'name': f'C{i+1:02d}', 'coord': cp, 'type': 'c3', 'charge': -0.12})

# Normalize total charge to -1
total_q = sum(a['charge'] for a in atoms)
for a in atoms: a['charge'] *= -1.0 / total_q
print(f"Total charge: {sum(a['charge'] for a in atoms):.3f}")

# Align distal carboxyl to MYR
fa_center = (o1d_pos + o2d_pos + ccx_pos) / 3
fa_dir = (c_pos[-1] - fa_center) / np.linalg.norm(c_pos[-1] - fa_center)

v = np.cross(fa_dir, chain_dir_u); s = np.linalg.norm(v); c = np.dot(fa_dir, chain_dir_u)
if s > 1e-8:
    vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    R = np.eye(3) + vx + vx@vx*((1-c)/(s*s))
else:
    R = np.eye(3)

for a in atoms:
    a['coord'] = (a['coord'] - fa_center) @ R.T + carboxyl_center

# Verify placement
for n in ('O1D','O2D'):
    a = next(x for x in atoms if x['name']==n)
    d = np.linalg.norm(a['coord'] - myr_o1)
    print(f"  {n} -> MYR O1: {d:.1f} A")

# Write mol2 (heavy atoms only)
n_heavy = len(atoms)
with open("c18_monoacid_fa3.mol2", 'w') as f:
    f.write("@<TRIPOS>MOLECULE\nFAH\n")
    f.write(f" {n_heavy} 0 0 0 0\nSMALL\nGAFF2\nC18 monoacid at FA3\n\n")
    f.write("@<TRIPOS>ATOM\n")
    for a in atoms:
        f.write(f"{a['aid']:6d} {a['name']:6s} {a['coord'][0]:10.4f} {a['coord'][1]:10.4f} {a['coord'][2]:10.4f} {a['type']:4s} 1 FAH {a['charge']:10.4f}\n")
    f.write("@<TRIPOS>BOND\n")
    # Chain C-C bonds
    bid = 0
    for i in range(n_c - 1):
        bid += 1; f.write(f"{bid:6d} {4+i:6d} {5+i:6d} 1\n")
    # Carboxyl bonds
    bid += 1; f.write(f"{bid:6d} {1:6d} {4:6d} 1\n")
    bid += 1; f.write(f"{bid:6d} {1:6d} {2:6d} 2\n")
    bid += 1; f.write(f"{bid:6d} {1:6d} {3:6d} 1\n")
    f.write("@<TRIPOS>SUBSTRUCTURE\n1 FAH 1 TEMP 0 **** **** 0 ROOT\n")

print(f"C18 monoacid mol2: {n_heavy} heavy atoms, ready for tleap")
