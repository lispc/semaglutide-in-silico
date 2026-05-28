#!/usr/bin/env python3
"""Build C18 diacid mol2 by extending the working monoacid template."""
import numpy as np

# ====== Build at origin (same as monoacid) ======
n_c = 18; bond_cc = 1.54; bond_ch = 1.09; ang_cc = 112*np.pi/180

c_pos = [np.array([0.,0.,0.]), np.array([bond_cc,0.,0.])]
for i in range(2, n_c):
    prev_u = (c_pos[-1] - c_pos[-2]) / np.linalg.norm(c_pos[-1] - c_pos[-2])
    new_d = np.cos(ang_cc)*prev_u + np.sin(ang_cc)*np.cross(np.array([0.,0.,1.]), prev_u)
    c_pos.append(c_pos[-1] + new_d * bond_cc)

# Perp helper
def perp_to(v):
    p = np.cross(v, np.array([0.,0.,1.]))
    if np.linalg.norm(p) < 0.01: p = np.cross(v, np.array([1.,0.,0.]))
    return p / np.linalg.norm(p)

# ====== Build ALL atoms (heavy + H) as a flat list ======
class Atom:
    def __init__(self, name, pos, atype, charge):
        self.name = name; self.pos = np.array(pos); self.atype = atype; self.charge = charge

heavy = []  # heavy atoms
hydro = []  # hydrogen atoms

# --- Distal carboxyl (in FA3 pocket) ---
to_c2 = (c_pos[1] - c_pos[0]) / np.linalg.norm(c_pos[1] - c_pos[0])
perp = perp_to(to_c2)
ccx_d = c_pos[0] - to_c2 * 1.52
heavy.append(Atom('C1D', ccx_d, 'c', 0.70))
heavy.append(Atom('O1D', ccx_d + perp * 1.08, 'o', -0.80))
heavy.append(Atom('O2D', ccx_d - perp * 1.08, 'o', -0.80))

# --- Carbon chain ---
for i, cp in enumerate(c_pos):
    heavy.append(Atom(f'C{i+1:02d}', cp, 'c3', -0.12))

# --- Proximal carboxyl (extends from C18) ---
to_c17 = (c_pos[-2] - c_pos[-1]) / np.linalg.norm(c_pos[-2] - c_pos[-1])
perp_p = perp_to(to_c17)
ccx_p = c_pos[-1] - to_c17 * 1.52
heavy.append(Atom('C1P', ccx_p, 'c', 0.70))
heavy.append(Atom('O1P', ccx_p + perp_p * 1.08, 'o', -0.80))
heavy.append(Atom('O2P', ccx_p - perp_p * 1.08, 'o', -0.80))

print(f"O1D-O2D: {np.linalg.norm(heavy[1].pos-heavy[2].pos):.2f}")
print(f"O1P-O2P: {np.linalg.norm(heavy[-2].pos-heavy[-1].pos):.2f}")

# --- Hydrogens on each chain carbon ---
# Heavy indices: 0=C1D, 1=O1D, 2=O2D, 3=C01, 4=C02, ..., 20=C18, 21=C1P, 22=O1P, 23=O2P
ci_start = 3  # index of C01 in heavy list
for i in range(n_c):
    cp = c_pos[i]
    if i == 0:
        v1 = heavy[0].pos - cp  # C1D
        v2 = c_pos[i+1] - cp
    elif i == n_c - 1:
        v1 = c_pos[i-1] - cp
        v2 = heavy[21].pos - cp  # C1P
    else:
        v1 = c_pos[i-1] - cp
        v2 = c_pos[i+1] - cp

    v1u = v1/np.linalg.norm(v1); v2u = v2/np.linalg.norm(v2)
    bis = -(v1u+v2u)
    if np.linalg.norm(bis) < 0.01: bis = np.array([1.,0.,0.])
    bis /= np.linalg.norm(bis)
    perp_h = np.cross(v1u, v2u)
    if np.linalg.norm(perp_h) < 0.01: perp_h = np.array([0.,1.,0.])
    perp_h /= np.linalg.norm(perp_h)
    hydro.append(Atom(f'H{i+1}A', cp + ((bis+perp_h)/np.linalg.norm(bis+perp_h))*bond_ch, 'hc', 0.06))
    hydro.append(Atom(f'H{i+1}B', cp + ((bis-perp_h)/np.linalg.norm(bis-perp_h))*bond_ch, 'hc', 0.06))

# Combine: heavy first, then hydrogens
all_atoms = heavy + hydro
n_heavy = len(heavy)
n_h = len(hydro)
n_total = len(all_atoms)

# Normalize charges to -2
total_q = sum(a.charge for a in all_atoms)
for a in all_atoms: a.charge *= -2.0 / total_q
print(f"Total charge: {sum(a.charge for a in all_atoms):.3f}")
print(f"Heavy={n_heavy}, H={n_h}, Total={n_total}")

# ====== Align distal carboxyl to FA3 (MYR 1003) ======
myr_o1 = np.array([13.4, 8.3, 10.2]); myr_o2 = np.array([14.2, 10.0, 9.0])
myr_c1 = np.array([13.9, 9.5, 10.1]); myr_c14 = np.array([9.6, 12.6, 17.6])
carb_c = (myr_o1+myr_o2+myr_c1)/3
chain_u = (myr_c14-myr_c1)/np.linalg.norm(myr_c14-myr_c1)

fa_center = (heavy[0].pos + heavy[1].pos + heavy[2].pos) / 3  # C1D + O1D + O2D
fa_dir = (c_pos[-1] - fa_center) / np.linalg.norm(c_pos[-1] - fa_center)

v = np.cross(fa_dir, chain_u); s = np.linalg.norm(v); c = np.dot(fa_dir, chain_u)
if s > 1e-8:
    vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    R = np.eye(3) + vx + vx@vx*((1-c)/(s*s))
else:
    R = np.eye(3)

for a in all_atoms:
    a.pos = (a.pos - fa_center) @ R.T + carb_c

# ====== Write mol2 ======
with open("c18_diacid_fa3.mol2", 'w') as f:
    f.write("@<TRIPOS>MOLECULE\nFAH\n")
    f.write(f" {n_total} 0 0 0 0\nSMALL\nGAFF2\nC18 diacid at FA3\n\n")

    f.write("@<TRIPOS>ATOM\n")
    for aid, a in enumerate(all_atoms):
        f.write(f"{aid+1:6d} {a.name:6s} {a.pos[0]:10.4f} {a.pos[1]:10.4f} {a.pos[2]:10.4f} {a.atype:4s} 1 FAH {a.charge:10.4f}\n")

    f.write("@<TRIPOS>BOND\n")
    bid = 0
    # C1D to O1D, O2D, C01 (heavy indices 0→1,2,3)
    bid+=1; f.write(f"{bid:6d} {1:6d} {2:6d} 2\n")  # C1D=O1D (double bond)
    bid+=1; f.write(f"{bid:6d} {1:6d} {3:6d} 1\n")  # C1D-O2D
    bid+=1; f.write(f"{bid:6d} {1:6d} {4:6d} 1\n")  # C1D-C01

    # C-C chain: heavy indices 3→4, 4→5, ..., 19→20 (C01→C02, ..., C17→C18)
    for i in range(n_c - 1):
        bid+=1; f.write(f"{bid:6d} {4+i:6d} {5+i:6d} 1\n")

    # C18 → C1P: heavy index 20→21
    bid+=1; f.write(f"{bid:6d} {21:6d} {22:6d} 1\n")  # C18 → C1P

    # C1P to O1P, O2P: heavy indices 21→22, 23
    bid+=1; f.write(f"{bid:6d} {22:6d} {23:6d} 2\n")  # C1P=O1P
    bid+=1; f.write(f"{bid:6d} {22:6d} {24:6d} 1\n")  # C1P-O2P

    # C-H bonds: hydrogens start at index n_heavy (0-based) = n_heavy + 1 (1-based)
    for i, a in enumerate(all_atoms):
        if a.atype != 'hc': continue
        hn = a.name  # e.g., "H01A"
        cnum = int(hn[1:].rstrip('AB'))  # extract carbon number
        carbon_1based = 4 + cnum - 1  # heavy index of this carbon (1-based)
        bid+=1; f.write(f"{bid:6d} {carbon_1based:6d} {i+1:6d} 1\n")

    f.write("@<TRIPOS>SUBSTRUCTURE\n1 FAH 1 TEMP 0 **** **** 0 ROOT\n")

print("Saved: c18_diacid_fa3.mol2")
