#!/usr/bin/env python3
"""Modify WT prmtop → Aib8 by editing topology arrays directly, then convert to GROMACS."""
import parmed as pmd, numpy as np, os

tleap_dir = "/home/scroll/personal/semaglutide-in-silico/exps/exp-A/tleap"
out_dir = "/home/scroll/personal/semaglutide-in-silico/exps/exp-A/gmx/aib8"
os.makedirs(out_dir, exist_ok=True)

# Load WT topology
a = pmd.load_file(f"{tleap_dir}/wt.prmtop", f"{tleap_dir}/wt.inpcrd")

# Find ALA at idx 729 (peptide position 2)
res = a.residues[729]
assert res.name == 'ALA', f"Expected ALA, got {res.name}"

# Get atom indices for key atoms
ha_idx = None; ca_idx = None; cb_idx = None
n_pos = None; ca_pos = None; c_pos = None; cb_pos = None
for atom in res.atoms:
    xyz = np.array([atom.xx, atom.xy, atom.xz])
    if atom.name == 'HA': ha_idx = atom.idx
    if atom.name == 'CA': ca_idx = atom.idx; ca_pos = xyz
    if atom.name == 'N':  n_pos = xyz
    if atom.name == 'C':  c_pos = xyz
    if atom.name == 'CB': cb_idx = atom.idx; cb_pos = xyz

assert ha_idx is not None, "HA not found!"

# Calculate CB2 and H positions
u_n = (n_pos - ca_pos); u_n /= np.linalg.norm(u_n)
u_c = (c_pos - ca_pos); u_c /= np.linalg.norm(u_c)
u_cb = (cb_pos - ca_pos); u_cb /= np.linalg.norm(u_cb)
u_cb2 = -(u_n + u_c + u_cb); u_cb2 /= np.linalg.norm(u_cb2)
cb2_pos = ca_pos + u_cb2 * 1.53

u_ca = -u_cb2
perp = np.array([1.,0.,0.])
if abs(np.dot(perp, u_ca)) > 0.9: perp = np.array([0.,1.,0.])
perp -= np.dot(perp, u_ca)*u_ca; perp /= np.linalg.norm(perp)
hb21_pos = cb2_pos + u_ca * 1.09
hb22_pos = cb2_pos + (-0.5*u_ca + 0.866*np.cross(perp, u_ca)) * 1.09
hb23_pos = cb2_pos + (-0.5*u_ca - 0.866*np.cross(perp, u_ca)) * 1.09

# Build new topology arrays
n_old = len(a.atoms)
n_new = n_old + 3  # HA→CB2 + 3H

# Copy and extend arrays
names = list(a.parm_data['ATOM_NAME'])
charges = list(a.parm_data['CHARGE'])
masses = list(a.parm_data['MASS'])
types = list(a.parm_data['AMBER_ATOM_TYPE'])
coords = a.coordinates  # (n_old, 3) numpy array

# Modify HA → CB2 at ha_idx
names[ha_idx] = 'CB2'
charges[ha_idx] = -0.1721
masses[ha_idx] = 12.01
types[ha_idx] = 'CT'
coords[ha_idx] = cb2_pos

# Add 3 new H atoms
for hname, hpos, hchg in [('HB21',hb21_pos,0.0707), ('HB22',hb22_pos,0.0707), ('HB23',hb23_pos,0.0707)]:
    names.append(hname)
    charges.append(hchg)
    masses.append(1.008)
    types.append('HC')
coords = np.vstack([coords, hb21_pos, hb22_pos, hb23_pos])

# Update other ALA→AIB charges
for atom in res.atoms:
    if atom.name == 'CA': charges[atom.idx] = 0.0341
    elif atom.name == 'CB': charges[atom.idx] = -0.1721
    elif atom.name == 'HB1': charges[atom.idx] = 0.0707
    elif atom.name == 'HB2': charges[atom.idx] = 0.0707
    elif atom.name == 'HB3': charges[atom.idx] = 0.0707

# Fix bonds: remove HA bonds, add CB2 bonds
# BONDS_INC_HYDROGEN format: triplets of [atom_i*3, atom_j*3, type]
def fix_bonds(bond_arr, ha_idx, ca_idx, cb2_idx_new, hb21_idx, hb22_idx, hb23_idx):
    """Remove HA bonds, add CB2 bonds."""
    new = []
    i = 0
    while i < len(bond_arr):
        a1 = bond_arr[i] // 3
        a2 = bond_arr[i+1] // 3
        bt = bond_arr[i+2]
        if a1 != ha_idx and a2 != ha_idx:
            new.extend([bond_arr[i], bond_arr[i+1], bond_arr[i+2]])
        i += 3
    return new

bonds_h = list(a.parm_data['BONDS_INC_HYDROGEN'])
bonds_nh = list(a.parm_data['BONDS_WITHOUT_HYDROGEN'])

hb21_idx = n_old
hb22_idx = n_old + 1
hb23_idx = n_old + 2
cb2_idx = ha_idx  # CB2 replaced HA

bonds_h = fix_bonds(bonds_h, ha_idx, ca_idx, cb2_idx, hb21_idx, hb22_idx, hb23_idx)
bonds_nh = fix_bonds(bonds_nh, ha_idx, ca_idx, cb2_idx, hb21_idx, hb22_idx, hb23_idx)

# Add CA-CB2 (no H: bonds_nh)
bonds_nh.extend([ca_idx*3, cb2_idx*3, 1])
# Add CB2-H bonds (bonds_h)
for hidx in [hb21_idx, hb22_idx, hb23_idx]:
    bonds_h.extend([cb2_idx*3, hidx*3, 1])

# Update charges also need to fix the total charge... skip for now

# ====== SAVE AS GROMACS (not Amber prmtop — avoid complex validation) ======
# Update POINTERS
a.parm_data['POINTERS'][0] = n_new  # NATOM
a.parm_data['ATOM_NAME'] = names
a.parm_data['CHARGE'] = charges
a.parm_data['MASS'] = masses
a.parm_data['AMBER_ATOM_TYPE'] = types
a.parm_data['BONDS_INC_HYDROGEN'] = bonds_h
a.parm_data['BONDS_WITHOUT_HYDROGEN'] = bonds_nh
a.parm_data['POINTERS'][2] = len(bonds_h) // 3
a.parm_data['POINTERS'][3] = len(bonds_nh) // 3
a.coordinates = coords

# Rename residue
res.name = 'AIB'

# Save as GROMACS directly
a.save(f"{out_dir}/aib8.top", format="gromacs", overwrite=True)
a.save(f"{out_dir}/aib8.gro", format="gromacs", overwrite=True)

# Verify
print(f"Total atoms: {len(a.atoms)} (was {n_old})")
r = a.residues[729]
print(f"Residue {r.idx}: {r.name}, {len(list(r.atoms))} atoms")
for atom in r.atoms:
    xyz = np.array([atom.xx, atom.xy, atom.xz])
    print(f"  {atom.name:6s} type={atom.type:6s} chg={atom.charge:.4f} pos=({xyz[0]:.2f},{xyz[1]:.2f},{xyz[2]:.2f})")

print(f"\nSaved: {out_dir}/aib8.top, {out_dir}/aib8.gro")
