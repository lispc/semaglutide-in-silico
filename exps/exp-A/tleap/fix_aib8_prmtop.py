#!/usr/bin/env python3
"""
Properly modify the WT prmtop to Aib8 by directly editing topology arrays.
This works at the Amber parameter file level, bypassing ParmEd's residue API.
"""
import parmed as pmd, numpy as np, copy, os

label = "aib8"
tleap_dir = "/home/scroll/personal/semaglutide-in-silico/exps/exp-A/tleap"
out_dir = f"/home/scroll/personal/semaglutide-in-silico/exps/exp-A/md/{label}"

# Load WT
a = pmd.load_file(f"{tleap_dir}/wt.prmtop", f"{tleap_dir}/wt.inpcrd")

# Find ALA at idx 729
res = a.residues[729]
assert res.name == 'ALA'

# Collect atom indices and current data
ala_atoms = list(res.atoms)
print(f"ALA idx 729: {len(ala_atoms)} atoms")

# Find HA atom
ha_idx = None; ha_atom = None
ca_atom = None; cb_atom = None
for atom in ala_atoms:
    if atom.name == 'HA': ha_idx = atom.idx; ha_atom = atom
    if atom.name == 'CA': ca_atom = atom
    if atom.name == 'CB': cb_atom = atom

if ha_atom is None:
    print("ERROR: HA not found")
    exit(1)

# Calculate CB2 position
n_pos = None; ca_pos = None; c_pos = None; cb_pos = None
for atom in res.atoms:
    xyz = np.array([atom.xx, atom.xy, atom.xz])
    if atom.name == 'N': n_pos = xyz
    if atom.name == 'CA': ca_pos = xyz
    if atom.name == 'C': c_pos = xyz
    if atom.name == 'CB': cb_pos = xyz

u_n = (n_pos - ca_pos); u_n /= np.linalg.norm(u_n)
u_c = (c_pos - ca_pos); u_c /= np.linalg.norm(u_c)
u_cb = (cb_pos - ca_pos); u_cb /= np.linalg.norm(u_cb)
u_cb2 = -(u_n + u_c + u_cb); u_cb2 /= np.linalg.norm(u_cb2)
cb2_pos = ca_pos + u_cb2 * 1.53

u_ca = -u_cb2
perp = np.array([1.0,0.0,0.0])
if abs(np.dot(perp, u_ca)) > 0.9: perp = np.array([0.0,1.0,0.0])
perp -= np.dot(perp, u_ca)*u_ca; perp /= np.linalg.norm(perp)
cos120, sin120 = -0.5, 0.866
hb21_pos = cb2_pos + u_ca * 1.09
hb22_pos = cb2_pos + (cos120*u_ca + sin120*np.cross(perp, u_ca))*1.09
hb23_pos = cb2_pos + (cos120*u_ca - sin120*np.cross(perp, u_ca))*1.09

# Build new topology arrays
# We need to:
# 1. Replace HA with CB2 (same index slot — reuse HA's position in arrays)
# 2. Add 3 new atoms (HB21, HB22, HB23) at the end of the topology
# 3. Update ATOM_NAME, CHARGE, MASS, ATOM_TYPE_INDEX, etc.
# 4. Update BONDS arrays
# 5. Add coordinates for new atoms

n_atoms_old = len(a.atoms)
n_atoms_new = n_atoms_old + 3  # HA → CB2 + 3 new H atoms

# Work with ParmEd's internal arrays
# Atom names
old_names = list(a.parm_data['ATOM_NAME'])
old_charges = list(a.parm_data['CHARGE'])
old_masses = list(a.parm_data['MASS'])
old_atom_types = list(a.parm_data['AMBER_ATOM_TYPE'])
old_coords = a.coordinates

# Replace HA → CB2 at ha_idx
old_names[ha_idx] = 'CB2'
old_charges[ha_idx] = -0.1721
old_masses[ha_idx] = 12.01
old_atom_types[ha_idx] = 'CT'
old_coords[ha_idx] = cb2_pos

# Add HB21, HB22, HB23
for hname, hpos in [('HB21', hb21_pos), ('HB22', hb22_pos), ('HB23', hb23_pos)]:
    old_names.append(hname)
    old_charges.append(0.0707)
    old_masses.append(1.008)
    old_atom_types.append('HC')
    old_coords = np.vstack([old_coords, hpos])

# Update CA charge
for i, atom in enumerate(res.atoms):
    if atom.name == 'CA':
        old_charges[atom.idx] = 0.0341
    elif atom.name == 'CB':
        old_charges[atom.idx] = -0.1721
    elif atom.name == 'HB1':
        old_charges[atom.idx] = 0.0707
    elif atom.name == 'HB2':
        old_charges[atom.idx] = 0.0707
    elif atom.name == 'HB3':
        old_charges[atom.idx] = 0.0707

# Update POINTERS first (before coordinates, which validates against NATOM)
a.parm_data['POINTERS'][0] = n_atoms_new  # NATOM

# Update ParmEd internal arrays
a.parm_data['ATOM_NAME'] = old_names
a.parm_data['CHARGE'] = old_charges
a.parm_data['MASS'] = old_masses
a.parm_data['AMBER_ATOM_TYPE'] = old_atom_types
a.coordinates = old_coords

# Update residue name
res.name = 'AIB'

# Update bonds: remove HA bonds, add CB2 bonds
# Remove bonds involving HA
bonds_inc = list(a.parm_data['BONDS_INC_HYDROGEN'])
bonds_no_h = list(a.parm_data['BONDS_WITHOUT_HYDROGEN'])

# Find and remove HA bonds
def remove_bonds_with_atom(bond_list, atom_idx):
    new_bonds = []
    i = 0
    while i < len(bond_list):
        a1 = bond_list[i] // 3
        a2 = bond_list[i+1] // 3
        idx = bond_list[i+2]  # bond type index
        if a1 != atom_idx and a2 != atom_idx:
            new_bonds.extend([bond_list[i], bond_list[i+1], bond_list[i+2]])
        i += 3
    return new_bonds

bonds_inc = remove_bonds_with_atom(bonds_inc, ha_idx)
bonds_no_h = remove_bonds_with_atom(bonds_no_h, ha_idx)

# Add CA-CB2 bond (heavy-heavy = BONDS_WITHOUT_HYDROGEN)
# For bonds with hydrogen, bond atom1 = atom_with_H, atom2 = H
# For bonds without hydrogen, bond atom1 = atom1, atom2 = atom2
# Each bond entry is [atom1*3, atom2*3, bond_type_idx]

# CA-CB2 (no hydrogen involved)
ca_idx = ca_atom.idx
cb2_idx = ha_idx  # CB2 replaced HA at this index
bonds_no_h.extend([ca_idx*3, cb2_idx*3, 1])  # bond_type 1 = single bond

# CB2-H bonds (hydrogen bonds)
hb21_idx = n_atoms_old      # first new atom
hb22_idx = n_atoms_old + 1
hb23_idx = n_atoms_old + 2
for h_idx in [hb21_idx, hb22_idx, hb23_idx]:
    bonds_inc.extend([cb2_idx*3, h_idx*3, 1])

a.parm_data['BONDS_INC_HYDROGEN'] = bonds_inc
a.parm_data['BONDS_WITHOUT_HYDROGEN'] = bonds_no_h

# Update bond counts in POINTERS
# POINTERS[2] = NBONH (bonds with hydrogen), POINTERS[3] = MBONA (bonds without H)
a.parm_data['POINTERS'][2] = len(bonds_inc) // 3
a.parm_data['POINTERS'][3] = len(bonds_no_h) // 3

# Save modified prmtop/inpcrd
os.makedirs(out_dir, exist_ok=True)
a.save(f"{out_dir}/aib8_mod2.prmtop", overwrite=True, format='amber')
a.save(f"{out_dir}/aib8_mod2.inpcrd", overwrite=True, format='rst7')

# Verify
a2 = pmd.load_file(f"{out_dir}/aib8_mod2.prmtop", f"{out_dir}/aib8_mod2.inpcrd")
r = a2.residues[729]
print(f"Residue {r.idx}: {r.name}, {len(list(r.atoms))} atoms")
for atom in r.atoms:
    xyz = np.array([atom.xx, atom.xy, atom.xz])
    print(f"  {atom.name:6s} type={atom.type:6s} chg={atom.charge:.4f} pos=({xyz[0]:.2f},{xyz[1]:.2f},{xyz[2]:.2f})")
print(f"\nTotal atoms: {len(a2.atoms)} (was {n_atoms_old})")
print(f"Saved: {out_dir}/aib8_mod2.prmtop")
