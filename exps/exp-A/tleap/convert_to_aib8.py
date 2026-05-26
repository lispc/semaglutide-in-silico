#!/usr/bin/env python3
"""Convert WT prmtop/inpcrd to Aib8 by modifying peptide residue 2 (ALA→AIB)."""
import parmed as pmd
import numpy as np
import shutil, os

# Load WT
wt_dir = "/home/scroll/personal/semaglutide-in-silico/exps/exp-A/tleap"
a = pmd.load_file(f"{wt_dir}/wt.prmtop", f"{wt_dir}/wt.inpcrd")

# Find peptide ALA (idx 729, number 729, the 2nd residue of GLP-1 chain P)
res = a.residues[729]
assert res.name == 'ALA', f"Expected ALA at idx 729, got {res.name}"

# Get coordinates
ca_pos = None; cb_pos = None; n_pos = None; c_pos = None; ha_atom = None
for atom in res.atoms:
    xyz = np.array([atom.xx, atom.xy, atom.xz])
    if atom.name == 'CA': ca_pos = xyz
    if atom.name == 'CB': cb_pos = xyz
    if atom.name == 'N': n_pos = xyz
    if atom.name == 'C': c_pos = xyz
    if atom.name == 'HA': ha_atom = atom

# Calculate CB2 position (tetrahedral, opposite to N+C+CB from CA)
u_n = (n_pos - ca_pos); u_n /= np.linalg.norm(u_n)
u_c = (c_pos - ca_pos); u_c /= np.linalg.norm(u_c)
u_cb = (cb_pos - ca_pos); u_cb /= np.linalg.norm(u_cb)
u_cb2 = -(u_n + u_c + u_cb); u_cb2 /= np.linalg.norm(u_cb2)
cb2_pos = ca_pos + u_cb2 * 1.53

# CB2 hydrogens
u_ca = -u_cb2
perp = np.array([1.0,0.0,0.0])
if abs(np.dot(perp, u_ca)) > 0.9: perp = np.array([0.0,1.0,0.0])
perp -= np.dot(perp, u_ca) * u_ca; perp /= np.linalg.norm(perp)
cos120, sin120 = -0.5, 0.866
hb21_pos = cb2_pos + u_ca * 1.09
hb22_pos = cb2_pos + (cos120*u_ca + sin120*np.cross(perp, u_ca)) * 1.09
hb23_pos = cb2_pos + (cos120*u_ca - sin120*np.cross(perp, u_ca)) * 1.09

# Check for clashes with nearby atoms (DPP-4 residues within 5A)
print(f"CB2 position: {cb2_pos}")
nearby = []
for atom in a.atoms:
    xyz = np.array([atom.xx, atom.xy, atom.xz])
    d = np.linalg.norm(cb2_pos - xyz)
    if d < 2.0 and atom.residue.idx != 729:
        nearby.append((d, atom.name, atom.residue.name, atom.residue.number))
nearby.sort()
if nearby:
    print(f"WARNING: {len(nearby)} atoms within 2A of CB2:")
    for d, an, rn, ri in nearby[:5]:
        print(f"  {d:.1f}A: {an}@{rn}{ri}")
else:
    print("No clashes detected for CB2")

# Modify topology
# 1. Remove HA
res.delete_atom(ha_atom)

# 2. Add CB2 atom (CT type, charge -0.1721)
cb2 = pmd.Atom(name='CB2', type='CT', charge=-0.1721, mass=12.01)
res.add_atom(cb2)

# 3. Add HB21, HB22, HB23 (HC type, charge 0.0707)
for hname in ['HB21', 'HB22', 'HB23']:
    h = pmd.Atom(name=hname, type='HC', charge=0.0707, mass=1.008)
    res.add_atom(h)

# 4. Update charges
charges = {
    'CA': 0.0341,    # was 0.0337
    'CB': -0.1721,   # was -0.1825
    'HB1': 0.0707,   # was 0.0603
    'HB2': 0.0707,   # was 0.0603
    'HB3': 0.0707,   # was 0.0603
}
for atom in res.atoms:
    if atom.name in charges:
        atom.charge = charges[atom.name]

# 5. Rename residue
res.name = 'AIB'

# 6. Add bonds for CB2
atom_map = {a.name: a for a in res.atoms}
# ParmEd uses bonds list, add Bond objects
import parmed as pmd_module
a.bonds.append(pmd_module.Bond(atom_map['CA'], atom_map['CB2']))
for hname in ['HB21', 'HB22', 'HB23']:
    a.bonds.append(pmd_module.Bond(atom_map['CB2'], atom_map[hname]))

print(f"Modified residue: {res.name} with {len(res.atoms)} atoms")
for atom in res.atoms:
    print(f"  {atom.name:6s} type={atom.type:6s} chg={atom.charge:.4f}")

# Update coordinates
# Find CB2 and its H atoms (last 4 atoms added)
new_atoms = list(res.atoms)
cb2_atom = new_atoms[-4]   # CB2
hb21_atom = new_atoms[-3]  # HB21
hb22_atom = new_atoms[-2]  # HB22
hb23_atom = new_atoms[-1]  # HB23

# Set coordinates (need to use atom.xx etc.)
cb2_atom.xx, cb2_atom.xy, cb2_atom.xz = cb2_pos[0], cb2_pos[1], cb2_pos[2]
hb21_atom.xx, hb21_atom.xy, hb21_atom.xz = hb21_pos[0], hb21_pos[1], hb21_pos[2]
hb22_atom.xx, hb22_atom.xy, hb22_atom.xz = hb22_pos[0], hb22_pos[1], hb22_pos[2]
hb23_atom.xx, hb23_atom.xy, hb23_atom.xz = hb23_pos[0], hb23_pos[1], hb23_pos[2]

# Save modified topology
out_dir = "/home/scroll/personal/semaglutide-in-silico/exps/exp-A/md/aib8"
os.makedirs(out_dir, exist_ok=True)
a.save(f"{out_dir}/aib8_modified.prmtop", overwrite=True)
a.save(f"{out_dir}/aib8_modified.inpcrd", overwrite=True)
print(f"\nSaved: {out_dir}/aib8_modified.prmtop, aib8_modified.inpcrd")
