#!/usr/bin/env python3
"""Fix carboxylate O-O distances in mol2 file."""
import numpy as np

with open("c18_diacid_fa3.mol2") as f:
    lines = f.readlines()

atom_start = None; atom_end = None
for i, line in enumerate(lines):
    if "@<TRIPOS>ATOM" in line: atom_start = i + 1
    if "@<TRIPOS>BOND" in line and atom_start is not None: atom_end = i; break

atoms = []
for i in range(atom_start, atom_end):
    if not lines[i].strip(): continue
    parts = lines[i].split()
    if len(parts) >= 6:
        atoms.append({'idx': i, 'aid': int(parts[0]), 'name': parts[1],
                      'coord': np.array([float(parts[2]), float(parts[3]), float(parts[4])]),
                      'type': parts[5]})

for n1, n2, c_name in [('O1P', 'O2P', 'C01'), ('O1D', 'O2D', 'C18X')]:
    a1 = next((a for a in atoms if a['name'] == n1), None)
    a2 = next((a for a in atoms if a['name'] == n2), None)
    c = next((a for a in atoms if a['name'] == c_name), None)
    if not all([a1, a2, c]): continue
    d = np.linalg.norm(a1['coord'] - a2['coord'])
    print(f"{n1}-{n2}: {d:.2f} A")
    if d < 1.5:
        mid = (a1['coord'] + a2['coord']) / 2
        to_c = c['coord'] - mid
        perp = np.cross(to_c, np.array([0., 0., 1.]))
        if np.linalg.norm(perp) < 0.01:
            perp = np.cross(to_c, np.array([1., 0., 0.]))
        perp = perp / np.linalg.norm(perp)
        a1['coord'] = mid + perp * 1.08
        a2['coord'] = mid - perp * 1.08
        nd = np.linalg.norm(a1['coord'] - a2['coord'])
        print(f"  Fixed: {nd:.2f} A")
        for a in [a1, a2]:
            parts = lines[a['idx']].split()
            lines[a['idx']] = f"{a['aid']:6d} {a['name']:6s} {a['coord'][0]:10.4f} {a['coord'][1]:10.4f} {a['coord'][2]:10.4f} {a['type']:4s} 1 FAH {float(parts[-1]):10.4f}\n"

with open("c18_diacid_fa3.mol2", 'w') as f:
    f.writelines(lines)
print("Fixed mol2 saved.")
