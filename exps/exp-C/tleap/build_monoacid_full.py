#!/usr/bin/env python3
"""Build C18 monoacid mol2 WITH hydrogens, positioned at FA3 site."""
import numpy as np

# ====== Build at origin ======
n_c = 18; bond_cc = 1.54; bond_ch = 1.09; ang_cc = 112*np.pi/180

c_pos = [np.array([0.,0.,0.]), np.array([bond_cc,0.,0.])]
for i in range(2, n_c):
    prev_u = (c_pos[-1] - c_pos[-2]) / np.linalg.norm(c_pos[-1] - c_pos[-2])
    new_d = np.cos(ang_cc)*prev_u + np.sin(ang_cc)*np.cross(np.array([0.,0.,1.]), prev_u)
    c_pos.append(c_pos[-1] + new_d * bond_cc)

# Carboxyl at C1
to_c2 = (c_pos[1] - c_pos[0]) / np.linalg.norm(c_pos[1] - c_pos[0])
perp = np.cross(to_c2, np.array([0.,0.,1.]))
perp = perp / (np.linalg.norm(perp) or 1)
ccx = c_pos[0] - to_c2 * 1.52
o1d = ccx + perp * 1.08; o2d = ccx - perp * 1.08

atoms = []; aid = 0
def add(name, pos, atype, charge):
    global aid; aid += 1
    atoms.append({'aid': aid, 'name': name, 'coord': pos, 'type': atype, 'charge': charge})
    return aid

# Carboxyl group
add('C1X', ccx, 'c', 0.70); add('O1D', o1d, 'o', -0.80); add('O2D', o2d, 'o', -0.80)

# Carbon chain
c_ids = []
for i, cp in enumerate(c_pos):
    c_ids.append(add(f'C{i+1:02d}', cp, 'c3', -0.12))

# H atoms on each carbon
for i in range(n_c):
    cp = c_pos[i]; ci = c_ids[i]
    if i == n_c - 1:  # CH3 terminal
        prev_u = (c_pos[i-1] - cp) / np.linalg.norm(c_pos[i-1] - cp)
        p1 = np.array([0.,0.,1.]); p2 = np.cross(prev_u, p1)
        if np.linalg.norm(p2) < 0.01: p2 = np.cross(prev_u, np.array([1.,0.,0.]))
        p2 /= np.linalg.norm(p2); p1 = np.cross(p2, prev_u)
        base = cp + prev_u * bond_ch * 0.333
        add(f'H{i+1}A', base + p2*bond_ch*0.943, 'hc', 0.06)
        add(f'H{i+1}B', base + (-0.5*p2+0.866*p1)*bond_ch*0.943, 'hc', 0.06)
        add(f'H{i+1}C', base + (-0.5*p2-0.866*p1)*bond_ch*0.943, 'hc', 0.06)
    else:  # CH2
        if i == 0: v1 = ccx - cp; v2 = c_pos[i+1] - cp
        else: v1 = c_pos[i-1] - cp; v2 = c_pos[i+1] - cp
        v1u = v1/np.linalg.norm(v1); v2u = v2/np.linalg.norm(v2)
        bis = -(v1u+v2u)
        if np.linalg.norm(bis) < 0.01: bis = np.array([1.,0.,0.])
        bis /= np.linalg.norm(bis)
        perp = np.cross(v1u, v2u)
        if np.linalg.norm(perp) < 0.01: perp = np.array([0.,1.,0.])
        perp /= np.linalg.norm(perp)
        add(f'H{i+1}A', cp + ((bis+perp)/np.linalg.norm(bis+perp))*bond_ch, 'hc', 0.06)
        add(f'H{i+1}B', cp + ((bis-perp)/np.linalg.norm(bis-perp))*bond_ch, 'hc', 0.06)

# Normalize charges
total_q = sum(a['charge'] for a in atoms)
for a in atoms: a['charge'] *= -1.0 / total_q

# ====== Align to FA3 ======
myr_o1 = np.array([13.4, 8.3, 10.2]); myr_o2 = np.array([14.2, 10.0, 9.0])
myr_c1 = np.array([13.9, 9.5, 10.1]); myr_c14 = np.array([9.6, 12.6, 17.6])
carb_c = (myr_o1+myr_o2+myr_c1)/3
chain_u = (myr_c14-myr_c1)/np.linalg.norm(myr_c14-myr_c1)

fa_center = (o1d+o2d+ccx)/3
fa_dir = (c_pos[-1]-fa_center)/np.linalg.norm(c_pos[-1]-fa_center)

v = np.cross(fa_dir, chain_u); s = np.linalg.norm(v); c = np.dot(fa_dir, chain_u)
R = np.eye(3) + (np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]) if s>1e-8 else 0) + \
    (np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])@np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]]))*((1-c)/(s*s)) if s>1e-8 else 0

for a in atoms:
    a['coord'] = (a['coord'] - fa_center) @ R.T + carb_c

# Write mol2
n_atoms = len(atoms)
with open("c18_monoacid_fa3.mol2", 'w') as f:
    f.write("@<TRIPOS>MOLECULE\nFAH\n")
    f.write(f" {n_atoms} 0 0 0 0\nSMALL\nGAFF2\nC18 monoacid at FA3 with H\n\n")
    f.write("@<TRIPOS>ATOM\n")
    for a in atoms:
        f.write(f"{a['aid']:6d} {a['name']:6s} {a['coord'][0]:10.4f} {a['coord'][1]:10.4f} {a['coord'][2]:10.4f} {a['type']:4s} 1 FAH {a['charge']:10.4f}\n")
    f.write("@<TRIPOS>BOND\n")
    bid = 0
    for i in range(n_c-1):
        bid+=1; f.write(f"{bid:6d} {4+i:6d} {5+i:6d} 1\n")  # C-C
    bid+=1; f.write(f"{bid:6d} {1:6d} {4:6d} 1\n")  # C1X-C01
    bid+=1; f.write(f"{bid:6d} {1:6d} {2:6d} 2\n")  # C1X=O1D
    bid+=1; f.write(f"{bid:6d} {1:6d} {3:6d} 1\n")  # C1X-O2D
    # C-H bonds
    for a in atoms:
        if a['name'].startswith('H'):
            hn = a['name']
            cnum = int(hn[1:].rstrip('ABC'))
            bid+=1; f.write(f"{bid:6d} {4+cnum-1:6d} {a['aid']:6d} 1\n")
    f.write("@<TRIPOS>SUBSTRUCTURE\n1 FAH 1 TEMP 0 **** **** 0 ROOT\n")

print(f"Built: {n_atoms} atoms (with H), total charge {sum(a['charge'] for a in atoms):.3f}")
print(f"Heavy atoms: {sum(1 for a in atoms if a['type']!='hc')}")
print(f"Hydrogens: {sum(1 for a in atoms if a['type']=='hc')}")
